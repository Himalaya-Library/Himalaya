// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "himalaya/mh2_eft/ThresholdCalculator.hpp"
#include "himalaya/mh2_eft/EFTFlags.hpp"
#include "himalaya/mh2_eft/threshold_loop_functions.hpp"
#include "himalaya/misc/Constants.hpp"
#include "himalaya/misc/Li2.hpp"
#include "himalaya/misc/Logger.hpp"
#include "himalaya/misc/Powers.hpp"
#include <cmath>
#include <complex>
#include <stdexcept>

/**
 * @file ThresholdCalculator.cpp
 *
 * @brief Implementation of threshold corrections class to express the
 * Higgs mass calculation in terms of SM MS-bar / MSSM DR'-bar
 * parameters.
 */

namespace himalaya {
namespace mh2_eft {

namespace {

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      return std::fabs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      return is_zero(a - b, prec);
   }

   template <typename T>
   bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
         return true;

      if (std::abs(a) < std::numeric_limits<T>::epsilon() ||
          std::abs(b) < std::numeric_limits<T>::epsilon())
         return false;

      return std::abs((a - b)/a) < prec;
   }

   double calc_cw(double mW, double mZ)
   {
      return std::abs(mZ) > std::numeric_limits<double>::epsilon()
         ? mW/mZ
         : 1.0;
   }

/// threshold loop functions from [1504.05200] Eq.(14)
double f1HD(double x) noexcept
{
   const double x2 = pow2(x);
   const double den = 1 - x2;

   if (std::abs(den) < 1e-5)
      return -1;

   return x2/den * std::log(x2);
}

/// threshold loop functions from [1504.05200] Eq.(15)
double f2HD(double x) noexcept
{
   const double x2 = pow2(x);
   const double den = 1 - x2;

   if (std::abs(den) < 1e-5)
      return 0.5;

   return (1 + x2/den * std::log(x2))/den;
}

/// threshold loop functions from [1504.05200] Eq.(17)
double f3HD(double x) noexcept
{
   using himalaya::dilog;

   const double x2 = pow2(x);
   const double den = 1 - x2;
   const std::complex<double> cden = den;

   if (std::abs(den) < 1e-5)
      return -9./4;

   const double logx2 = std::log(x2);

   const std::complex<double> result =
      (-1 + x2*(2 + 2*x2))/pow2(den) * (
         logx2*std::log(cden) + dilog(x2) - z2 - x2*logx2);

   return result.real();
}


bool isfinite(double exact, double shifted, double limit) noexcept
{
   // checks if the threshold correction in the general mass case is
   // finite and smaller than 1
   if (!std::isfinite(exact) || !std::isfinite(shifted)
       || (std::abs(exact) > 1. && std::abs(shifted) > 1.))
      return false;

   // checks if the difference of the shifted result to the limit is
   // greater than the difference of the exact result to the limit
   // which should indicate the divergence of the general mass case
   if (std::abs(shifted - limit) >= std::abs(exact - limit))
      return false;

   return true;
}

double deltaxyz(double x, double y, double z) noexcept
{
   return threshold_loop_functions::delta_xyz(x, y, z);
}

/**
 * \f$\Phi(x,y,z)\f$ function.  The arguments x, y and z are
 * interpreted as squared masses. Taken from FlexibleSUSY.
 *
 * Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
 *
 * @param x squared mass
 * @param y squared mass
 * @param z squared mass
 *
 * @return \f$\Phi(x,y,z)\f$
 */
double phixyz(double x, double y, double z) noexcept
{
   return threshold_loop_functions::phi_xyz(x, y, z);
}

} // anonymous namespace

using namespace threshold_loop_functions;

/**
 * Constructor
 * @param p_ a HimalayaInterface struct
 * @param msq2_ the averaged squark mass of the first two generations squared
 * @param verbose a bool enable the output of the parameter validation. Enabled by default
 * @param check a boolean which indicates if the threshold corrections should be tested
 */
ThresholdCalculator::ThresholdCalculator(
   const Parameters& p_, double msq2_, bool verbose, bool check)
   : p(p_), msq2(msq2_)
{
   p.validate(verbose);

   if (!std::isfinite(msq2_))
      msq2 = p.calculateMsq2();

   if (!check) {
      // Set mass limit for threshold corrections
      const double mQ3 = std::sqrt(p.mq2(2,2));
      const double mU3 = std::sqrt(p.mu2(2,2));
      const double m3 = p.MG;
      const double eps = mQ3*0.01;
      const double eps2 = mU3*0.01;
      const double msq2Save = msq2;
      const double v = std::sqrt(pow2(p.vu) + pow2(p.vd));
      const double pref = sqrt2*p.Mt*pow4(p.g3/(4*Pi))/v;

      if (std::abs(mQ3-mU3) < eps && std::abs(mU3 - m3) < eps && std::abs(m3 - std::sqrt(msq2)) < eps) {
         const double lim = pref*getDeltaYtAlphas2(Limits::DEGENERATE, 1);
         const double exact = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         p.mu2(2,2) = pow2(mQ3 + std::abs(mU3 - mQ3)/2.);
         p.MG = mQ3 + std::abs(m3 - mQ3)/2.;
         msq2 = pow2(mQ3 + std::abs(std::sqrt(msq2) - mQ3)/2.);
         const double exactShifted = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         if (!isfinite(exact, exactShifted, lim)) p.massLimit3LThreshold = static_cast<int>(Limits::DEGENERATE);
      } else if (std::abs(mQ3 - mU3) < eps && std::abs(mU3 - m3) < eps) {
         const double lim = pref*getDeltaYtAlphas2(Limits::MQ3_EQ_MU3_EQ_M3, 1);
         const double exact = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         p.mu2(2,2) = pow2(mQ3 + std::abs(mU3 - mQ3)/2.);
         p.MG = mQ3 + std::abs(m3 - mQ3)/2.;
         const double exactShifted = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         if (!isfinite(exact, exactShifted, lim)) p.massLimit3LThreshold = static_cast<int>(Limits::MQ3_EQ_MU3_EQ_M3);
      } else if (std::abs(mQ3 - mU3) < eps) {
         const double lim = pref*getDeltaYtAlphas2(Limits::MQ3_EQ_MU3, 1);
         const double exact = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         p.mu2(2,2) = pow2(mQ3 + std::abs(mU3 - mQ3)/2.);
         p.MG = mQ3 + std::abs(m3 - mQ3)/2.;
         const double exactShifted = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         if (!isfinite(exact, exactShifted, lim)) p.massLimit3LThreshold = static_cast<int>(Limits::MQ3_EQ_MU3);
      } else if (std::abs(mQ3 - m3) < eps) {
         const double lim = pref*getDeltaYtAlphas2(Limits::MQ3_EQ_M3, 1);
         const double exact = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         p.mu2(2,2) = pow2(mQ3 + std::abs(mU3 - mQ3)/2.);
         p.MG = mQ3 + std::abs(m3 - mQ3)/2.;
         const double exactShifted = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         if (!isfinite(exact, exactShifted, lim)) p.massLimit3LThreshold = static_cast<int>(Limits::MQ3_EQ_M3);
      } else if (std::abs(mU3 - m3) < eps2) {
         const double lim = pref*getDeltaYtAlphas2(Limits::MU3_EQ_M3, 1);
         const double exact = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         p.mu2(2,2) = pow2(mQ3 + std::abs(mU3 - mQ3)/2.);
         p.MG = mQ3 + std::abs(m3 - mQ3)/2.;
         const double exactShifted = pref*getDeltaYtAlphas2(Limits::GENERAL, 1);
         if (!isfinite(exact, exactShifted, lim)) p.massLimit3LThreshold = static_cast<int>(Limits::MU3_EQ_M3);
      }

      // reset possible parameter shifts
      p.mu2(2,2) = pow2(mU3);
      p.MG = m3;
      msq2 = msq2Save;
   }
}

/**
 * Returns a specific threshold corrections for a given mass limit
 * @param variable coupling order of threshold correction
 * @param scheme renormalization scheme. Choices are {DRbar'} only!
 * @param omitLogs an integer to omit all log mu terms
 * @return a threshold correction for a given variable in a given scheme for a suitable mass limit
 */
double ThresholdCalculator::getThresholdCorrection(
   ThresholdCouplingOrders variable, RenSchemes scheme, int omitLogs) const
{
   double thresholdCorrection = 0.;
   const auto limit = static_cast<Limits>(p.massLimit3LThreshold);

   if (scheme != RenSchemes::TEST && scheme != RenSchemes::DRBARPRIME) {
      INFO_MSG("Your renormalization scheme is not compatible with the"
               " implemented threshold corrections!");
   }

   switch (variable) {
      case ThresholdCouplingOrders::G3_AS:{
         thresholdCorrection = getDeltaG3Alphas(omitLogs);
         switch (scheme) {
            case RenSchemes::DRBARPRIME:
               thresholdCorrection = - thresholdCorrection;
               break;
            default:
               break;
         }
      }
      break;
      case ThresholdCouplingOrders::YT_AS:{
         thresholdCorrection = getDeltaYtAlphas(limit, omitLogs);
         switch (scheme) {
            case RenSchemes::DRBARPRIME:
               thresholdCorrection = - thresholdCorrection;
               break;
            default:
               break;
         }
      }
      break;
      case ThresholdCouplingOrders::YT_AS2:{
         thresholdCorrection = getDeltaYtAlphas2(limit, omitLogs);
         switch (scheme) {
            case RenSchemes::DRBARPRIME:{
               const double dytas = getDeltaYtAlphas(limit, omitLogs);
               thresholdCorrection = 2 * getDeltaG3Alphas(omitLogs) * dytas
                  + pow2(dytas) - thresholdCorrection;
            }
            break;
         default:
            break;
         }
      }
      break;
      case ThresholdCouplingOrders::LAMBDA_AT:{
         thresholdCorrection = getDeltaLambdaAlphat(limit, omitLogs);
      }
      break;
      case ThresholdCouplingOrders::LAMBDA_AT_AS:{
         thresholdCorrection = getDeltaLambdaAlphatAlphas(limit, omitLogs);
         switch (scheme) {
            case RenSchemes::DRBARPRIME:
               thresholdCorrection = thresholdCorrection
                  - 4 * getDeltaLambdaAlphat(limit, omitLogs) *
                  getDeltaYtAlphas(limit, omitLogs);
               break;
            default:
               break;
         }
      }
      break;
      // Note that the genuine contribution of lambda_atas2 is unknown and thus
      // set to 0 (note: here are the reconstructed DR' logs included)
      // The lines below just convert it from MSbar to DRbar
      case ThresholdCouplingOrders::LAMBDA_AT_AS2:{
         thresholdCorrection = getDeltaLambdaAlphatAlphas2(limit, omitLogs);
         switch (scheme) {
            case RenSchemes::DRBARPRIME:{
               const double dg3as = getDeltaG3Alphas(omitLogs);
               const double dytas = getDeltaYtAlphas(limit, omitLogs);
               thresholdCorrection = thresholdCorrection +
               (getDeltaLambdaAlphatAlphas(limit, omitLogs)
               * (-2 * dg3as - 4 * dytas)
               + getDeltaLambdaAlphat(limit, omitLogs) * (8 * dg3as * dytas
                  + 10 * pow2(dytas) - 4 * getDeltaYtAlphas2(limit, omitLogs)));
            }
            break;
         default:
            break;
         }
      }
      break;
      case ThresholdCouplingOrders::LAMBDA_YB2_G12:{
         return getDeltaLambdaYb2G12(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_G14:{
         return getDeltaLambdaG14(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_REG_G14:{
         return getDeltaLambdaRegG14();
      }
      case ThresholdCouplingOrders::LAMBDA_CHI_G14:{
         return getDeltaLambdaChiG14(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_CHI_G24:{
         return getDeltaLambdaChiG24(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_G24:{
         return getDeltaLambdaG24(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_REG_G24:{
         return getDeltaLambdaRegG24();
      }
      case ThresholdCouplingOrders::LAMBDA_G12_G22:{
         return getDeltaLambdaG12G22(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_REG_G12_G22:{
         return getDeltaLambdaRegG12G22();
      }
      case ThresholdCouplingOrders::LAMBDA_CHI_G12_G22:{
         return getDeltaLambdaChiG12G22(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YB2_G22:{
         return getDeltaLambdaYb2G22(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YB4:{
         return getDeltaLambdaYb4(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YT2_G12:{
         return getDeltaLambdaYt2G12(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YT2_G22:{
         return getDeltaLambdaYt2G22(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YTAU2_G12:{
         return getDeltaLambdaYtau2G12(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YTAU2_G22:{
         return getDeltaLambdaYtau2G22(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YTAU4:{
         return getDeltaLambdaYtau4(omitLogs);
      }
      case ThresholdCouplingOrders::G1_G1:{
         return getDeltaG1G1(omitLogs);
      }
      case ThresholdCouplingOrders::G2_G2:{
         return getDeltaG2G2(omitLogs);
      }
      case ThresholdCouplingOrders::VEV_YT2:{
         return getDeltaVevYt2(limit);
      }
      case ThresholdCouplingOrders::YT_YB:{
         return getDeltaYtYb(omitLogs);
      }
      case ThresholdCouplingOrders::YT_YT:{
         return getDeltaYtYt(omitLogs);
      }
      case ThresholdCouplingOrders::YTAU_YTAU:{
         return getDeltaYtauYtau(omitLogs);
      }
      case ThresholdCouplingOrders::YB_YB:{
        return getDeltaYbYb(omitLogs);
      }
      case ThresholdCouplingOrders::YB_YT:{
        return getDeltaYbYt(omitLogs);
      }
      case ThresholdCouplingOrders::YB_AS:{
        return getDeltaYbAs(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YB4_G32:{
         return getDeltaLambdaYb4G32(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YB6:{
         return getDeltaLambdaYb6(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YT6:{
         return getDeltaLambdaYt6(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YTAU6:{
         return getDeltaLambdaYtau6(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YT2_YB4:{
         return getDeltaLambdaYt2Yb4(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YT4_YB2:{
         return getDeltaLambdaYt4Yb2(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YTAU4_YB2:{
        return getDeltaLambdaYtau4Yb2(omitLogs);
      }
      case ThresholdCouplingOrders::LAMBDA_YTAU2_YB4:{
        return getDeltaLambdaYtau2Yb4(omitLogs);
      }
      case ThresholdCouplingOrders::VEV_G12:{
         return getDeltaVevG12(omitLogs);
      }
      case ThresholdCouplingOrders::VEV_G22:{
         return getDeltaVevG22(omitLogs);
      }
      case ThresholdCouplingOrders::VEV_YB2:{
         return getDeltaVevYb2();
      }
      case ThresholdCouplingOrders::VEV_YTAU2:{
         return getDeltaVevYtau2();
      }
      case ThresholdCouplingOrders::YTAU_YB:{
         return getDeltaYtauYb();
      }
      default:
         break;
   };

   return thresholdCorrection;
}

double ThresholdCalculator::getDeltaLambdaYb2G12(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mD3 = std::sqrt(mD32);
   const double Xb2 = pow2(p.Ad(2,2) - p.mu*p.vu/p.vd);
   const double c2beta = std::cos(2*atan(p.vu/p.vd));
   const double lmD3MR = omitLogs*log(mD32 / MR2);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);

   return -3*c2beta*(Xb2*(3*(F3(mQ3/mD3) + F4(mQ3/mD3)) + 2*c2beta*F5(mQ3/mD3)) +
      4*mD3*mQ3*(2*lmD3MR + lmQ3MR))/(20.*mD3*mQ3);
}

double ThresholdCalculator::getDeltaLambdaG14(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);
   const double c4beta = std::cos(4*beta);
   const double c8beta = std::cos(8*beta);
   const double s4beta2 = pow2(sin(4*beta));
   const double lmAMR = omitLogs*log(pow2(p.MA) / MR2);
   const double lmD1MR = omitLogs*log(p.md2(0,0) / MR2);
   const double lmD2MR = omitLogs*log(p.md2(1,1) / MR2);
   const double lmD3MR = omitLogs*log(p.md2(2,2) / MR2);
   const double lmE1MR = omitLogs*log(p.me2(0,0) / MR2);
   const double lmE2MR = omitLogs*log(p.me2(1,1) / MR2);
   const double lmE3MR = omitLogs*log(p.me2(2,2) / MR2);
   const double lmL1MR = omitLogs*log(p.ml2(0,0) / MR2);
   const double lmL2MR = omitLogs*log(p.ml2(1,1) / MR2);
   const double lmL3MR = omitLogs*log(p.ml2(2,2) / MR2);
   const double lmQ1MR = omitLogs*log(p.mq2(0,0) / MR2);
   const double lmQ2MR = omitLogs*log(p.mq2(1,1) / MR2);
   const double lmQ3MR = omitLogs*log(p.mq2(2,2) / MR2);
   const double lmU1MR = omitLogs*log(p.mu2(0,0) / MR2);
   const double lmU2MR = omitLogs*log(p.mu2(1,1) / MR2);
   const double lmU3MR = omitLogs*log(p.mu2(2,2) / MR2);

   return (-108*s4beta2 - 3*(-29 + 4*c4beta + 9*c8beta)*lmAMR + 16*pow2(c2beta)*
      (2*(lmD1MR + lmD2MR + lmD3MR + 3*(lmE1MR + lmE2MR + lmE3MR)) + 3*(lmL1MR +
      lmL2MR + lmL3MR) + lmQ1MR + lmQ2MR + lmQ3MR + 8*(lmU1MR + lmU2MR + lmU3MR)
      ))/800.;
}

double ThresholdCalculator::getDeltaLambdaRegG14() const
{
   return -9/50.;
}

double ThresholdCalculator::getDeltaLambdaChiG14(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double c4be = std::cos(4*beta);
   const double s2be = std::sin(2*beta);
   const double s6be = std::sin(6*beta);
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);

   return (-3*(22*lmUMR - 2*c4be*lmUMR + 7*(3 + c4be)*f1(p.M1/p.mu) - 9*(-1
      + c4be)*f3(p.M1/p.mu) + 31*s2be*F5(p.M1/p.mu) - 2*F7(p.M1/p.mu)
      - 2*c4be*F7(p.M1/p.mu) - F5(p.M1/p.mu)*s6be))/200.;
}

double ThresholdCalculator::getDeltaLambdaChiG24(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double c4beta = std::cos(4*beta);
   const double s2beta = std::sin(2*beta);
   const double s6beta = std::sin(6*beta);
   const double lm2MR = omitLogs*log(pow2(p.M2) / MR2);
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);

   return (-82*lmUMR - 10*c4beta*lmUMR - 27*(3 + c4beta)*f2(p.M2/p.mu) 
        - 93*s2beta*F5(p.M2/p.mu)
        + 6*F7(p.M2/p.mu) + 6*c4beta*F7(p.M2/p.mu) - 8*lm2MR - 8*c4beta*lm2MR
        - 42*f4(p.M2/p.mu)*pow2(s2beta) + 3*F5(p.M2/p.mu)*s6beta)/24.;
}

double ThresholdCalculator::getDeltaLambdaG24(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);
   const double c4beta = std::cos(4*beta);
   const double c8beta = std::cos(8*beta);
   const double s4beta = std::sin(4*beta);
   const double lmAMR = omitLogs*log(pow2(p.MA) / MR2);
   const double lmL1MR = omitLogs*log(p.ml2(0,0) / MR2);
   const double lmL2MR = omitLogs*log(p.ml2(1,1) / MR2);
   const double lmL3MR = omitLogs*log(p.ml2(2,2) / MR2);
   const double lmQ1MR = omitLogs*log(p.mq2(0,0) / MR2);
   const double lmQ2MR = omitLogs*log(p.mq2(1,1) / MR2);
   const double lmQ3MR = omitLogs*log(p.mq2(2,2) / MR2);

   return (-36*pow2(s4beta) + (53 - 28*c4beta - 9*c8beta)*lmAMR + 16*pow2(c2beta)
      *(lmL1MR + lmL2MR + lmL3MR + 3*(lmQ1MR + lmQ2MR + lmQ3MR)))/96.;
}

double ThresholdCalculator::getDeltaLambdaRegG24() const
{
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);

   return (-9 + 2*pow2(c2beta))/6.;
}

double ThresholdCalculator::getDeltaLambdaG12G22(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double c4beta = std::cos(4*beta);
   const double c8beta = std::cos(8*beta);
   const double s4beta = std::sin(4*beta);
   const double lmAMR = omitLogs*log(pow2(p.MA) / MR2);

   return -3*(12*pow2(s4beta) + (-7 + 4*c4beta + 3*c8beta)*lmAMR)/80.;
}

double ThresholdCalculator::getDeltaLambdaRegG12G22() const
{
   return -3/5.;
}

double ThresholdCalculator::getDeltaLambdaChiG12G22(int omitLogs) const
{
   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);
   const double c4beta = std::cos(4*beta);
   const double s2beta = std::sin(2*beta);
   const double s6beta = std::sin(6*beta);
   const double m1mu = p.M1/p.mu;
   const double m2mu = p.M2/p.mu;
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);

   return (-24*lmUMR + 24*c4beta*lmUMR + 3*s2beta*F5(m2mu) - 42*f6(m1mu,m2mu) - 14*
        c4beta*f6(m1mu,m2mu) + 6*F7(m2mu) + 6*c4beta*F7(m2mu) - 64*s2beta*f8(m1mu
        ,m2mu) + 4*F7(m1mu)*pow2(c2beta) - 32*f5(m1mu, m2mu)*pow2(s2beta) -
        4*f7(m1mu,m2mu)*pow2(s2beta) + 3*F5(m2mu)*s6beta + F5(m1mu)*(
        s2beta + s6beta))/40.;
}

double ThresholdCalculator::getDeltaLambdaYb2G22(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double Xb2 = pow2(p.Ad(2,2) - p.mu*p.vu/p.vd);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mD3 = std::sqrt(mD32);
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);

   return -(c2beta*Xb2*(3*F4(mQ3/mD3) + c2beta*F5(mQ3/mD3)))/(2.*mD3*mQ3)
      - 3*c2beta*lmQ3MR;
}

double ThresholdCalculator::getDeltaLambdaYb4(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double Xb2 = pow2(p.Ad(2,2) - p.mu*p.vu/p.vd);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mD3 = std::sqrt(mD32);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);
   const double lmD3MR = omitLogs*log(mD32 / MR2);

   return (12*mD3*mQ3*Xb2*F1(mQ3/mD3) - pow2(Xb2)*F2(mQ3/mD3))/(mD32*mQ32)
      + 6*(lmD3MR + lmQ3MR);
}

double ThresholdCalculator::getDeltaLambdaYt2G12(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double Xt2 = pow2(p.Au(2,2) - p.mu*p.vd/p.vu);
   const double mQ32 = p.mq2(2,2);
   const double mU32 = p.mu2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mU3 = std::sqrt(mU32);
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);
   const double lmU3MR = omitLogs*log(mU32 / MR2);

   return (-3*c2beta*(-3*Xt2*F3(mQ3/mU3) + c2beta*Xt2*F5(mQ3/mU3)
      + 2*mQ3*mU3*(lmQ3MR - 4*lmU3MR)))/(10.*mQ3*mU3);
}

double ThresholdCalculator::getDeltaLambdaYt2G22(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double Xt2 = pow2(p.Au(2,2) - p.mu*p.vd/p.vu);
   const double mQ32 = p.mq2(2,2);
   const double mU32 = p.mu2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mU3 = std::sqrt(mU32);
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);

   return -c2beta*Xt2*(-3*F4(mQ3/mU3) + c2beta*F5(mQ3/mU3))/(2.*mQ3*mU3) + 3*
   c2beta*lmQ3MR;
}

double ThresholdCalculator::getDeltaLambdaYtau2G12(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double Xtau2 = pow2(p.Ae(2,2) - p.mu*p.vu/p.vd);
   const double mL32 = p.ml2(2,2);
   const double mE32 = p.me2(2,2);
   const double mL3 = std::sqrt(mL32);
   const double mE3 = std::sqrt(mE32);
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);
   const double lmE3MR = omitLogs*log(mE32 / MR2);
   const double lmL3MR = omitLogs*log(mL32 / MR2);

   return c2beta*((Xtau2*(-9*F3(mL3/mE3) + 3*F4(mL3/mE3) - 2*c2beta*F5(mL3/mE3))
   )/(mE3*mL3) + 12*(-2*lmE3MR + lmL3MR))/20.;
}

double ThresholdCalculator::getDeltaLambdaYtau2G22(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double Xtau2 = pow2(p.Ae(2,2) - p.mu*p.vu/p.vd);
   const double mL32 = p.ml2(2,2);
   const double mE32 = p.me2(2,2);
   const double mL3 = std::sqrt(mL32);
   const double mE3 = std::sqrt(mE32);
   const double beta = std::atan(p.vu/p.vd);
   const double c2beta = std::cos(2*beta);
   const double lmL3MR = omitLogs*log(mL32 / MR2);

   return -c2beta*Xtau2*(3*F4(mL3/mE3) + c2beta*F5(mL3/mE3))/(6.*mE3*mL3)
   - c2beta*lmL3MR;
}

double ThresholdCalculator::getDeltaLambdaYtau4(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double Xtau2 = pow2(p.Ae(2,2) - p.mu*p.vu/p.vd);
   const double mL32 = p.ml2(2,2);
   const double mE32 = p.me2(2,2);
   const double mL3 = std::sqrt(mL32);
   const double mE3 = std::sqrt(mE32);
   const double lmE3MR = omitLogs*log(mE32 / MR2);
   const double lmL3MR = omitLogs*log(mL32 / MR2);

   return (12*mE3*mL3*Xtau2*F1(mL3/mE3) - pow2(Xtau2)*F2(mL3/mE3))/(3.*mE32*mL32)
   + 2*(lmE3MR + lmL3MR);
}

double ThresholdCalculator::getDeltaG1G1(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double lmAMR = omitLogs*log(pow2(p.MA) / MR2);
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);
   const double lmD1MR = omitLogs*log(p.md2(0,0) / MR2);
   const double lmD2MR = omitLogs*log(p.md2(1,1) / MR2);
   const double lmD3MR = omitLogs*log(p.md2(2,2) / MR2);
   const double lmE1MR = omitLogs*log(p.me2(0,0) / MR2);
   const double lmE2MR = omitLogs*log(p.me2(1,1) / MR2);
   const double lmE3MR = omitLogs*log(p.me2(2,2) / MR2);
   const double lmL1MR = omitLogs*log(p.ml2(0,0) / MR2);
   const double lmL2MR = omitLogs*log(p.ml2(1,1) / MR2);
   const double lmL3MR = omitLogs*log(p.ml2(2,2) / MR2);
   const double lmQ1MR = omitLogs*log(p.mq2(0,0) / MR2);
   const double lmQ2MR = omitLogs*log(p.mq2(1,1) / MR2);
   const double lmQ3MR = omitLogs*log(p.mq2(2,2) / MR2);
   const double lmU1MR = omitLogs*log(p.mu2(0,0) / MR2);
   const double lmU2MR = omitLogs*log(p.mu2(1,1) / MR2);
   const double lmU3MR = omitLogs*log(p.mu2(2,2) / MR2);

   return (-3*lmAMR - 2*(lmD1MR + lmD2MR + lmD3MR + 3*(lmE1MR + lmE2MR + lmE3MR)
      ) - 3*(lmL1MR + lmL2MR + lmL3MR) - lmQ1MR - lmQ2MR - lmQ3MR - 12*lmUMR - 8
      *(lmU1MR + lmU2MR + lmU3MR))/60.;
}

double ThresholdCalculator::getDeltaG2G2(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double lm2MR = omitLogs*log(pow2(p.M2) / MR2);
   const double lmAMR = omitLogs*log(pow2(p.MA) / MR2);
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);
   const double lmL1MR = omitLogs*log(p.ml2(0,0) / MR2);
   const double lmL2MR = omitLogs*log(p.ml2(1,1) / MR2);
   const double lmL3MR = omitLogs*log(p.ml2(2,2) / MR2);
   const double lmQ1MR = omitLogs*log(p.mq2(0,0) / MR2);
   const double lmQ2MR = omitLogs*log(p.mq2(1,1) / MR2);
   const double lmQ3MR = omitLogs*log(p.mq2(2,2) / MR2);

   return (4 - 8*lm2MR - lmAMR - lmL1MR - lmL2MR - lmL3MR - 3*(lmQ1MR + lmQ2MR
      + lmQ3MR) - 4*lmUMR)/12.;
}

double ThresholdCalculator::getDeltaYtYt(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double cbeta = std::cos(beta);
   const double sbeta = std::sin(beta);
   const double Xt2 = pow2(p.Au(2,2) - p.mu*p.vd/p.vu);
   const double mQ3 = std::sqrt(p.mq2(2,2));
   const double mU3 = std::sqrt(p.mu2(2,2));
   const double lmAMR = omitLogs*log(pow2(p.MA) / MR2);
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);

   return (-2*pow2(sbeta)*Xt2*F5(mQ3/mU3)/(mQ3*mU3) + 8*F6(mQ3/p.mu) + 4
      *F6(mU3/p.mu) + pow2(cbeta)*(-3 + 6*lmAMR) + 6*lmUMR)/8.;
}

double ThresholdCalculator::getDeltaYtauYtau(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double cbeta = std::cos(beta);
   const double sbeta = std::sin(beta);
   const double Xtau2 = pow2(p.Ae(2,2) - p.mu*p.vu/p.vd);
   const double mL3 = std::sqrt(p.ml2(2,2));
   const double mE3 = std::sqrt(p.me2(2,2));
   const double lmAMR = omitLogs*log(pow2(p.MA) / MR2);
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);

   return - pow2(cbeta)*Xtau2*F5(mL3/mE3)/(12.*mE3*mL3) + (4*F6(mE3/p.mu) + 8
      *F6(mL3/p.mu) + pow2(sbeta)*(-3 + 6*lmAMR) + 6*lmUMR)/8.;
}

double ThresholdCalculator::getDeltaYbYb(int omitLogs) const
{
   const double beta = std::atan(p.vu / p.vd);
   const double cbeta = std::cos(beta);
   const double sbeta = std::sin(beta);
   const double mA2 = pow2(p.MA);
   const double MR2 = pow2(p.scale);
   const double mQ3 = std::sqrt(p.mq2(2,2));
   const double mD3 = std::sqrt(p.md2(2,2));
   const double Xb2 = pow2(p.Ad(2,2) - p.mu*p.vu/p.vd);
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);

   return (6*lmUMR + 4*F6(mD3/p.mu) + 8*F6(mQ3/p.mu) - (2*Xb2*F5(mQ3/mD3)*pow2(cbeta))/
        (mD3*mQ3*pow2(sbeta)) + 3*(-1 + 2*omitLogs*log(mA2/MR2))*pow2(sbeta))/8.;
}

double ThresholdCalculator::getDeltaYbYt(int omitLogs) const
{
   const double tbe = p.vu / p.vd;
   const double beta = std::atan(tbe);
   const double cbeta = std::cos(beta);
   const double sbeta = std::sin(beta);
   const double mA2 = pow2(p.MA);
   const double MR2 = pow2(p.scale);
   const double mU3 = std::sqrt(p.mu2(2,2));
   const double mQ3 = std::sqrt(p.mq2(2,2));
   const double mD3 = std::sqrt(p.md2(2,2));
   const double Xb2 = pow2(p.Ad(2,2) - p.mu*p.vu/p.vd);
   const double Xt = p.Au(2,2) - p.mu*p.vd/p.vu;
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);

   return (2*lmUMR + 4*F6(mU3/p.mu) + (4*tbe*Xt*F9(mQ3/p.mu,mU3/p.mu))/p.mu + (-1 + 2*omitLogs*log(
        mA2/MR2))*pow2(cbeta) - (2*Xb2*F5(mQ3/mD3)*pow2(cbeta))/(mD3*mQ3*pow2(
        sbeta)) - (4*pow2(Xt)*F5(mQ3/mU3)*pow2(sbeta))/(mQ3*mU3) + 8*(-1 + omitLogs*log(mA2/
        MR2))*pow2(sbeta))/8.;
}

double ThresholdCalculator::getDeltaYbAs(int omitLogs) const
{
   const double MR2 = pow2(p.scale);
   const double mD3 = std::sqrt(p.md2(2,2));
   const double mQ3 = std::sqrt(p.mq2(2,2));
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double m3 = p.MG;

   return (4*(1 + F6(mD3/m3) + F6(mQ3/m3) - (Xb*F9(mQ3/m3,mD3/m3))/m3 + omitLogs
        *log(pow2(m3)/MR2)))/3.;
}

double ThresholdCalculator::getDeltaYtYb(int omitLogs) const
{
   using std::log;

   const double MR2 = pow2(p.scale);
   const double beta = std::atan(p.vu/p.vd);
   const double cbeta = std::cos(beta);
   const double sbeta = std::sin(beta);
   const double mD3 = std::sqrt(p.md2(2,2));
   const double mQ3 = std::sqrt(p.mq2(2,2));
   const double mU3 = std::sqrt(p.mu2(2,2));
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double Xt = p.Au(2,2) - p.mu*p.vd/p.vu;
   const double lmAMR = omitLogs*log(pow2(p.MA) / MR2);
   const double lmUMR = omitLogs*log(pow2(p.mu) / MR2);

   return (-2*p.mu*mU3*pow2(Xb)*F5(mQ3/mD3)*pow2(cbeta)*(1 + pow2(sbeta)) + mD3*pow2(
        sbeta)*(-2*p.mu*pow2(Xt)*F5(mQ3/mU3)*pow2(sbeta) + mQ3*mU3*(4*p.mu*F6(mD3/p.mu) +
        4*p.vd/p.vu*Xb*F9(mQ3/p.mu,mD3/p.mu) + p.mu*(2*lmUMR - 8*pow2(cbeta) - pow2(
        sbeta) + 2*lmAMR*(4*pow2(cbeta) + pow2(sbeta))))))/(8.*mD3*mQ3*
        p.mu*mU3*pow2(sbeta));
}

double ThresholdCalculator::getDeltaVevYt2(Limits limit) const
{
   const double Xt2 = pow2(p.Au(2,2) - p.mu*p.vd/p.vu);
   const double mQ32 = p.mq2(2,2);
   const double mU32 = p.mu2(2,2);

   if (limit == Limits::MQ3_EQ_M3 || limit == Limits::MU3_EQ_M3) limit = Limits::GENERAL;

   switch (limit) {
      case Limits::GENERAL:
         return (3*Xt2*(-2*mQ32*mU32*log(mQ32/mU32) + pow2(mQ32) - pow2(mU32)))/(4.*pow3(
            mQ32 - mU32));
      default:
         return Xt2/(4.*mQ32);
   }
}

double ThresholdCalculator::getDeltaVevYb2() const
{
   const double Xb2 = pow2(p.Ad(2,2) - p.mu*p.vu/p.vd);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double eps = mQ32*0.01;

   const double exact = (3*Xb2*(2*mD32*mQ32*log(mQ32/mD32) + pow2(mD32)
      - pow2(mQ32)))/(4.*pow3(mD32 - mQ32));

   if (std::abs(mQ32 - mD32) < eps) {
     return Xb2/(4.*mQ32);
      /*const double mQ3 = std::sqrt(mQ32);
      const double mD3 = std::sqrt(mD32);
      const double lim = Xb2/(4.*mQ32);
      const double mD32Shifted2 = pow2(mQ3 + std::abs(mD3 - mQ3)/2.);
      const double shifted = (3*Xb2*(2*mD32Shifted2*mQ32*log(mQ32/mD32Shifted2)
        + pow2(mD32Shifted2) - pow2(mQ32)))/(4.*pow3(mD32Shifted2 - mQ32));
      std::cout << exact << " " << shifted << " " << lim << " party\n";
      if (!isfinite(exact, shifted, lim)) return lim;*/
   }

   return exact;
}

double ThresholdCalculator::getDeltaVevYtau2() const{
   const double Xtau2 = pow2(p.Ae(2,2) - p.mu*p.vu/p.vd);
   const double mL32 = p.ml2(2,2);
   const double mE32 = p.me2(2,2);
   const double eps = mL32*0.01;

   const double exact = (Xtau2*(2*log(mL32/mE32)*mE32*mL32 + pow2(mE32)
      - pow2(mL32)))/(4.*pow3(mE32 - mL32));

   if (std::abs(mL32 - mE32) < eps) {
     return Xtau2/(12.*mL32);
      /*const double mL3 = std::sqrt(mL32);
      const double mE3 = std::sqrt(mE32);
      const double lim = Xtau2/(12.*mL32);
      const double mE32Shifted2 = pow2(mL3 + std::abs(mE3 - mL3)/2.);
      const double shifted = (Xtau2*(2*log(mL32/mE32Shifted2)*mE32Shifted2*mL32
	 + pow2(mE32Shifted2) - pow2(mL32)))/(4.*pow3(mE32Shifted2 - mL32));
      if (!isfinite(exact, shifted, lim)) return lim;*/
   }

   return exact;
}

/**
 * @warning The implemented is currently incorrect, so don't trust the
 * result!
 */
double ThresholdCalculator::getDeltaVevG12(int omitLogs) const
{
   using std::acos;
   using std::log;
   using std::tan;

   const double sw = std::sin(acos(calc_cw(p.MW, p.MZ)));
   const double MR2 = pow2(p.scale);
   const double Mu = p.mu;
   const double Mu2 = pow2(p.mu);
   const double M1 = p.M1;
   const double mL3 = std::sqrt(p.ml2(2,2));
   const double mE3 = std::sqrt(p.me2(2,2));
   const double mL2 = std::sqrt(p.ml2(1,1));
   const double mE2 = std::sqrt(p.me2(1,1));
   const double mL1 = std::sqrt(p.ml2(0,0));
   const double mE1 = std::sqrt(p.me2(0,0));
   const double mQ3 = std::sqrt(p.mq2(2,2));
   const double mU3 = std::sqrt(p.mu2(2,2));
   const double mQ2 = std::sqrt(p.mq2(1,1));
   const double mU2 = std::sqrt(p.mu2(1,1));
   const double mQ1 = std::sqrt(p.mq2(0,0));
   const double mU1 = std::sqrt(p.mu2(0,0));
   const double mD3 = std::sqrt(p.mu2(2,2));
   const double mD2 = std::sqrt(p.mu2(1,1));
   const double mD1 = std::sqrt(p.mu2(0,0));
   const double beta = std::atan(p.vu/p.vd);
   const double sbeta = std::sin(beta);
   const double cbeta = std::cos(beta);
   const double s2beta = std::sin(2*beta);
   const double s6beta = std::sin(6*beta);
   const double c2beta = std::cos(2*beta);
   const double c4beta = std::cos(4*beta);
   const double eps = M1*0.01;
   const double lmUMR = omitLogs*log(pow2(Mu) / MR2);
   const double lmD3MR = omitLogs*log(pow2(mD3) / MR2);
   const double lmE3MR = omitLogs*log(pow2(mE3) / MR2);
   const double lmL3MR = omitLogs*log(pow2(mL3) / MR2);
   const double lmQ3MR = omitLogs*log(pow2(mQ3) / MR2);

   WARNING_MSG("Corrections to delta_v of O(g1^2) are currently incorrect!");

   if (std::abs(beta - Pi/4.) < Pi/4.*0.01)
      return 0.;

   if (std::abs(M1 - Mu) < eps) {
      const double lim = (pow2(sw)*(4096*lmD3MR + 12288*lmE3MR + 6144*lmL3MR + 2048*lmQ3MR +
        24576*lmUMR + 6144*log(pow2(p.MA)/MR2) + 16384*log(pow2(mU3)/MR2) + 4096*log(pow2(
        mD1)/MR2) + 4096*log(pow2(mD2)/MR2) + 12288*log(pow2(mE1)/MR2) + 12288*
        log(pow2(mE2)/MR2) + 6144*log(pow2(mL1)/MR2) + 6144*log(pow2(mL2)/MR2)
        + 2048*log(pow2(mQ1)/MR2) + 2048*log(pow2(mQ2)/MR2) + 16384*log(pow2(
        mU1)/MR2) + 16384*log(pow2(mU2)/MR2) - (80*(371 - 65*c4beta + 140*
        s2beta + 64*s6beta))/pow2(c2beta) + (15*(5985 - 2203*c4beta + 8196*
        s2beta)*log(MR2/pow2(p.M1)))/pow4(cbeta + sbeta)))/409600.;

      return lim;
   }

   const double exact =
       (pow2(sw)*(4*lmD3MR + 12*lmE3MR + 6*lmL3MR + 2*lmQ3MR + 24*lmUMR + 6*log(
        pow2(p.MA)/MR2) + 16*log(pow2(mU3)/MR2) + 4*log(pow2(mD1)/MR2) + 4*log(pow2(mD2)/
        MR2) + 12*log(pow2(mE1)/MR2) + 12*log(pow2(mE2)/MR2) + 6*log(pow2(mL1)/
        MR2) + 6*log(pow2(mL2)/MR2) + 2*log(pow2(mQ1)/MR2) + 2*log(pow2(mQ2)/
        MR2) + 16*log(pow2(mU1)/MR2) + 16*log(pow2(mU2)/MR2) + (30*log(MR2/Mu2)
        *pow3(Mu)*(6*M1*Mu2*s2beta + 2*(4 + c4beta)*Mu*pow2(M1) + (3*s2beta +
        s6beta)*pow3(M1) - 2*c4beta*pow3(Mu)))/(pow2(c2beta)*pow3(M1 - Mu)*
        pow3(M1 + Mu)) + (15*log(MR2/pow2(M1))*pow3(M1)*(-3*(5 + c4beta)*M1*Mu2
        - 12*Mu*s2beta*pow2(M1) + (-1 + 3*c4beta)*pow3(M1) - 2*(3*s2beta +
        s6beta)*pow3(Mu)))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) + (15*(-
        2*(9 + c4beta)*Mu2*pow2(M1) - Mu*(13*s2beta + s6beta)*pow3(M1) - M1*(
        13*s2beta + s6beta)*pow3(Mu) + 2*(-1 + c4beta)*pow4(M1) + (-1 + 3*
        c4beta)*pow4(Mu)))/(pow2(c2beta)*pow2(M1 - Mu)*pow2(M1 + Mu))))/400.;

   return exact;
}

/**
 * @warning The implemented is currently incorrect, so don't trust the
 * result!
 */
double ThresholdCalculator::getDeltaVevG22(int omitLogs) const
{
   using std::acos;
   using std::log;
   using std::tan;

   const double MR2 = pow2(p.scale);
   const double M1 = p.M1;
   const double M2 = p.M2;
   const double Mu = p.mu;
   const double Mu2 = pow2(Mu);
   const double cw = calc_cw(p.MW, p.MZ);
   const double sw = std::sin(acos(cw));
   const double beta = std::atan(p.vu/p.vd);
   const double cbeta = std::cos(beta);
   const double c2beta = std::cos(2*beta);
   const double c4beta = std::cos(4*beta);
   const double sbeta = std::sin(beta);
   const double s2beta = std::sin(2*beta);
   const double s6beta = std::sin(6*beta);
   const double eps = Mu*0.01;
   const double mL3 = std::sqrt(p.ml2(2,2));
   const double mL2 = std::sqrt(p.ml2(1,1));
   const double mL1 = std::sqrt(p.ml2(0,0));
   const double mQ3 = std::sqrt(p.mq2(2,2));
   const double mQ2 = std::sqrt(p.mq2(1,1));
   const double mQ1 = std::sqrt(p.mq2(0,0));
   const double lmUMR = omitLogs*log(pow2(Mu) / MR2);
   const double lmL3MR = omitLogs*log(pow2(mL3) / MR2);
   const double lmQ3MR = omitLogs*log(pow2(mQ3) / MR2);

   WARNING_MSG("Corrections to delta_v of O(g2^2) are currently incorrect!");

   if (std::abs(beta - Pi/4.) < Pi/4.*0.01)
      return 0.;

   if (std::abs(M1 - Mu) < eps && std::abs(M2 - Mu) < eps) {
      return (lmL3MR*pow2(cw))/24. + (lmQ3MR*pow2(cw))/8. + (lmUMR*pow2(cw))/6. + (
        log(pow2(p.MA)/MR2)*pow2(cw))/24. + (log(pow2(M2)/MR2)*pow2(cw))/3. + (log(
        pow2(mL1)/MR2)*pow2(cw))/24. + (log(pow2(mL2)/MR2)*pow2(cw))/24. + (
        log(pow2(mQ1)/MR2)*pow2(cw))/8. + (log(pow2(mQ2)/MR2)*pow2(cw))/8. - (
        998 + 408*s2beta + 256*s6beta + (101 - 428*s2beta - 64*s6beta)*pow2(cw)
        + 3*c4beta*(42 + 83*pow2(cw)))/(3072.*pow2(c2beta)) + (log(MR2/pow2(M1))
        *(24258 + 32776*s2beta - 3*(2123 + 2732*s2beta)*pow2(cw) + c4beta*(-
        8502 + 1819*pow2(cw))))/(16384.*pow4(cbeta + sbeta));
   }

   if (std::abs(M1 - Mu) < eps) {
      const double t2beta = tan(2*beta);
      const double c2w = cos(2*acos(cw));
      const double lim =
        (512*lmL3MR*pow2(cw) + 1536*lmQ3MR*pow2(cw) + 2048*lmUMR*pow2(cw) + 512*
        log(pow2(p.MA)/MR2)*pow2(cw) + 4096*log(pow2(M2)/MR2)*pow2(cw) + 512*log(pow2(
        mL1)/MR2)*pow2(cw) + 512*log(pow2(mL2)/MR2)*pow2(cw) + 1536*log(pow2(
        mQ1)/MR2)*pow2(cw) + 1536*log(pow2(mQ2)/MR2)*pow2(cw) + (2*log(MR2/
        pow2(M2))*pow3(M2)*(-3*(-2*Mu2*(605 + 385*c4beta + 36*s2beta)*pow3(M2)
        + 2*(-317 + 31*c4beta - 420*s2beta)*pow2(M2)*pow3(Mu) + Mu*(61 - 31*
        c4beta - 92*s2beta)*pow4(M2) + M2*(2013 + 2561*c4beta - 988*s2beta)*
        pow4(Mu) + (-35 + c4beta + 36*s2beta)*pow5(M2) + (-195 + 225*c4beta +
        676*s2beta + 768*s6beta)*pow5(Mu)) + pow2(cw)*((40*Mu2*pow2(cw)*pow3(M2
        + Mu)*pow4(cbeta + sbeta))/pow2(sw) + 3*(2*Mu2*(-985 + 67*c4beta - 100*
        s2beta)*pow3(M2) + 2*(-461 + 79*c4beta - 2916*s2beta)*pow2(M2)*pow3(Mu)
        + Mu*(49 - 43*c4beta - 92*s2beta)*pow4(M2) + M2*(-3755 + 1785*c4beta -
        1372*s2beta)*pow4(Mu) + (-35 + c4beta + 36*s2beta)*pow5(M2) + (-279 +
        269*c4beta - 1756*s2beta)*pow5(Mu)))))/(Mu2*pow2(c2beta)*pow3(M2 - Mu)*
        pow3(M2 + Mu)) - (log(MR2/Mu2)*(3*(-2*(974 + (-149 + 197*c2w)*c4beta +
        768*s6beta + 1496*c2beta*t2beta)*pow3(M2)*pow3(Mu) + Mu2*(3118 + (679 -
        863*c2w)*c4beta + 1496*c2beta*t2beta)*pow4(M2) + 2*(-3358 + 3*(-741 +
        365*c2w)*c4beta + 40*c2beta*t2beta)*pow2(M2)*pow4(Mu) + Mu*(165*(-1 +
        c2w)*c4beta + 2*(751 + 620*c2beta*t2beta))*pow5(M2) + M2*(3518 + (-613
        + 773*c2w)*c4beta + 4312*c2beta*t2beta)*pow5(Mu) + 256*(-2 + (-1 + c2w)
        *c4beta)*pow6(M2) + (4110 + 2007*c4beta + 497*c2w*c4beta + 2520*c2beta*
        t2beta)*pow6(Mu)) + 2*pow2(cw)*((20*pow2(cw)*pow3(Mu)*pow3(M2 + Mu)*(3*
        c2beta*s2beta*t2beta + 8*sbeta*pow3(cbeta) + 8*cbeta*pow3(sbeta) + 2*
        pow4(cbeta) + 2*pow4(sbeta)))/pow2(sw) - 3*(2*(-623 + 244*c2beta*
        t2beta)*pow3(M2)*pow3(Mu) + Mu2*(1475 + 572*c2beta*t2beta)*pow4(M2) +
        2*(1237 + 68*c2beta*t2beta)*pow2(M2)*pow4(Mu) + Mu*(751 + 620*c2beta*
        t2beta)*pow5(M2) + 3*M2*(709 + 2404*c2beta*t2beta)*pow5(Mu) - 256*pow6(
        M2) + (2547 + 1468*c2beta*t2beta)*pow6(Mu)))))/(pow2(c2beta)*pow3(M2 -
        Mu)*pow3(M2 + Mu)) + (2*((-3023 + 109*c4beta - 964*s2beta - 1024*
        s6beta)*pow3(Mu)*pow4(M2) + (-5947 - 7591*c4beta + 1300*s2beta + 1408*
        s6beta)*pow3(M2)*pow4(Mu) + 2*Mu2*(1255 + 403*c4beta - 148*s2beta - 64*
        s6beta)*pow5(M2) + 2*(4877 + 2321*c4beta + 1348*s2beta - 704*s6beta)*
        pow2(M2)*pow5(Mu) + 96*Mu*(-3 + c4beta + 4*s2beta)*pow6(M2) + 4*M2*(-
        319 + 545*c4beta - 416*s2beta + 256*s6beta)*pow6(Mu) - 3*(-35 + c4beta
        + 36*s2beta)*pow7(M2) + (-4139 + 529*c4beta - 4420*s2beta + 128*s6beta)
        *pow7(Mu) + pow2(cw)*((5555 - 1417*c4beta - 128*s6beta - 14780*c2beta*
        t2beta + 1024*pow2(c2beta))*pow3(Mu)*pow4(M2) + (-14645 + 5479*c4beta -
        256*s6beta + 13628*c2beta*t2beta + 2048*pow2(c2beta))*pow3(M2)*pow4(Mu)
        + (20*Mu2*pow2(cw)*pow2(M2 + Mu)*(pow3(M2) - pow3(Mu))*(3*pow2(c2beta)*
        pow2(t2beta) + 8*sbeta*pow3(cbeta) + 8*cbeta*pow3(sbeta) + 2*pow4(
        cbeta) + 2*pow4(sbeta)))/pow2(sw) - 2*Mu2*(2953 - 395*c4beta - 64*
        s6beta + 44*c2beta*t2beta + 512*pow2(c2beta))*pow5(M2) - 2*(-5533 +
        1343*c4beta - 128*s6beta + 8620*c2beta*t2beta + 1024*pow2(c2beta))*
        pow2(M2)*pow5(Mu) - 12*Mu*(-21 + 11*c4beta + 32*c2beta*t2beta)*pow6(M2)
        - 16*M2*(56 - 25*c4beta - 8*s6beta - 1031*c2beta*t2beta + 64*pow2(
        c2beta))*pow6(Mu) + 3*(-35 + c4beta + 36*c2beta*t2beta)*pow7(M2) + (
        6119 - 2917*c4beta - 128*s6beta + 4180*c2beta*t2beta + 1024*pow2(
        c2beta))*pow7(Mu))))/(Mu2*pow2(c2beta)*pow2(M2 + Mu)*pow3(M2 - Mu)))/
        12288.;

      return lim;
   }

   if (std::abs(M2 - Mu) < eps) {
      const double mA2 = pow2(p.MA);
      const double t2beta = tan(2*beta);
      const double w = acos(cw);
      const double lim = (-67584/pow2(c2beta) + 4096*lmL3MR*pow2(cw) + 12288*lmQ3MR*pow2(cw) +
        16384*lmUMR*pow2(cw) + 4096*log(mA2/MR2)*pow2(cw) + 32768*log(pow2(M2)/
        MR2)*pow2(cw) + 4096*log(pow2(mL1)/MR2)*pow2(cw) + 4096*log(pow2(mL2)/
        MR2)*pow2(cw) + 12288*log(pow2(mQ1)/MR2)*pow2(cw) + 12288*log(pow2(mQ2)
        /MR2)*pow2(cw) + (1464*Mu2*t2beta*pow2(1/cos(w))*pow3(M1))/(c2beta*pow2(
        M1 + Mu)*pow3(M1 - Mu)) - (17178*Mu2*pow2(1/cos(w))*pow3(M1))/(pow2(
        c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) - (32658*c4beta*Mu2*pow2(1/cos(w))*
        pow3(M1))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) - (15741*Mu2*cos(
        4*beta - 2*w)*pow2(1/cos(w))*pow3(M1))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(
        M1 - Mu)) - (15741*Mu2*cos(2*(2*beta + w))*pow2(1/cos(w))*pow3(M1))/(
        pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (21186*Mu2*(-1 + pow2(tan(
        w)))*pow3(M1))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) - (137640*
        t2beta*pow2(M1)*pow2(1/cos(w))*pow3(Mu))/(c2beta*pow2(M1 + Mu)*pow3(M1 -
        Mu)) - (20658*pow2(M1)*pow2(1/cos(w))*pow3(Mu))/(pow2(c2beta)*pow2(M1 +
        Mu)*pow3(M1 - Mu)) + (10614*c4beta*pow2(M1)*pow2(1/cos(w))*pow3(Mu))/(
        pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (5799*cos(4*beta - 2*w)*
        pow2(M1)*pow2(1/cos(w))*pow3(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 -
        Mu)) + (5799*cos(2*(2*beta + w))*pow2(M1)*pow2(1/cos(w))*pow3(Mu))/(pow2(
        c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (24090*pow2(M1)*(-1 + pow2(tan(
        w)))*pow3(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (393216*
        t2beta*pow3(M1)*pow3(Mu))/(c2beta*pow3(M1 - Mu)*pow3(M1 + Mu)) + (
        196608*pow3(M1)*pow3(Mu))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) +
        (131424*Mu*t2beta*pow2(1/cos(w))*pow4(M1))/(c2beta*pow2(M1 + Mu)*pow3(M1
        - Mu)) + (4392*Mu*pow2(1/cos(w))*pow4(M1))/(pow2(c2beta)*pow2(M1 + Mu)*
        pow3(M1 - Mu)) - (3192*c4beta*Mu*pow2(1/cos(w))*pow4(M1))/(pow2(c2beta)*
        pow2(M1 + Mu)*pow3(M1 - Mu)) - (1500*Mu*cos(4*beta - 2*w)*pow2(1/cos(w))*
        pow4(M1))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (3816*Mu*cos(2*
        w)*pow2(1/cos(w))*pow4(M1))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) -
        (1500*Mu*cos(2*(2*beta + w))*pow2(1/cos(w))*pow4(M1))/(pow2(c2beta)*pow2(
        M1 + Mu)*pow3(M1 - Mu)) - (233472*Mu2*pow4(M1))/(pow2(c2beta)*pow3(M1 -
        Mu)*pow3(M1 + Mu)) + (67584*c4beta*Mu2*pow4(M1))/(pow2(c2beta)*pow3(M1
        - Mu)*pow3(M1 + Mu)) + (18432*Mu2*s6beta*pow4(M1))/(pow2(c2beta)*pow3(
        M1 - Mu)*pow3(M1 + Mu)) - (552960*cbeta*Mu2*sbeta*pow4(M1))/(pow2(
        c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) - (20184*M1*t2beta*pow2(1/cos(w))*
        pow4(Mu))/(c2beta*pow2(M1 + Mu)*pow3(M1 - Mu)) - (42222*M1*pow2(1/cos(w))
        *pow4(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (33258*c4beta*
        M1*pow2(1/cos(w))*pow4(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) +
        (15993*M1*cos(4*beta - 2*w)*pow2(1/cos(w))*pow4(Mu))/(pow2(c2beta)*pow2(
        M1 + Mu)*pow3(M1 - Mu)) + (15993*M1*cos(2*(2*beta + w))*pow2(1/cos(w))*
        pow4(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (37926*M1*(-1 +
        pow2(tan(w)))*pow4(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (
        405504*pow2(M1)*pow4(Mu))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) -
        (92160*c4beta*pow2(M1)*pow4(Mu))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 +
        Mu)) - (18432*s6beta*pow2(M1)*pow4(Mu))/(pow2(c2beta)*pow3(M1 - Mu)*
        pow3(M1 + Mu)) + (749568*cbeta*sbeta*pow2(M1)*pow4(Mu))/(pow2(c2beta)*
        pow3(M1 - Mu)*pow3(M1 + Mu)) - (5760*t2beta*pow2(1/cos(w))*pow5(M1))/(
        c2beta*pow2(M1 + Mu)*pow3(M1 - Mu)) + (57984*pow2(1/cos(w))*pow5(M1))/(
        pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) - (384*c4beta*pow2(1/cos(w))*
        pow5(M1))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) - (144*cos(4*beta
        - 2*w)*pow2(1/cos(w))*pow5(M1))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)
        ) + (57696*cos(2*w)*pow2(1/cos(w))*pow5(M1))/(pow2(c2beta)*pow2(M1 + Mu)*
        pow3(M1 - Mu)) - (144*cos(2*(2*beta + w))*pow2(1/cos(w))*pow5(M1))/(pow2(
        c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) - (172032*Mu*pow5(M1))/(pow2(
        c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) + (24576*c4beta*Mu*pow5(M1))/(
        pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) - (6144*Mu*s6beta*pow5(M1))/(
        pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) - (602112*cbeta*Mu*sbeta*
        pow5(M1))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) + (2376*t2beta*
        pow2(1/cos(w))*pow5(Mu))/(c2beta*pow2(M1 + Mu)*pow3(M1 - Mu)) - (1515*
        cos(4*beta - 2*w)*pow2(1/cos(w))*pow5(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*
        pow3(M1 - Mu)) - (1515*cos(2*(2*beta + w))*pow2(1/cos(w))*pow5(Mu))/(
        pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (438*pow2(1/cos(w))*pow5(Mu))
        /(pow2(c2beta)*pow2(M1 + Mu)*pow3(-M1 + Mu)) + (1854*c4beta*pow2(1/cos(w)
        )*pow5(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(-M1 + Mu)) + (3570*(-1 +
        pow2(tan(w)))*pow5(Mu))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(-M1 + Mu)) - (
        24576*M1*pow5(Mu))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) - (24576*
        c4beta*M1*pow5(Mu))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) + (6144*
        M1*s6beta*pow5(Mu))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) - (
        184320*cbeta*M1*sbeta*pow5(Mu))/(pow2(c2beta)*pow3(M1 - Mu)*pow3(M1 +
        Mu)) + (768*c4beta*pow2(1/cos(w))*pow6(M1))/(Mu*pow2(c2beta)*pow2(M1 +
        Mu)*pow3(M1 - Mu)) + (6144*cbeta*sbeta*pow2(1/cos(w))*pow6(M1))/(Mu*pow2(
        c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (384*cos(4*beta - 2*w)*pow2(1/cos(
        w))*pow6(M1))/(Mu*pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (384*cos(
        2*(2*beta + w))*pow2(1/cos(w))*pow6(M1))/(Mu*pow2(c2beta)*pow2(M1 + Mu)*
        pow3(M1 - Mu)) + (2304*(-1 + pow2(tan(w)))*pow6(M1))/(Mu*pow2(c2beta)*
        pow2(M1 + Mu)*pow3(M1 - Mu)) + (2304*pow2(1/cos(w))*pow6(M1))/(Mu*pow2(
        c2beta)*pow2(M1 + Mu)*pow3(-M1 + Mu)) - (36864*pow6(M1))/(pow2(c2beta)*
        pow3(M1 - Mu)*pow3(M1 + Mu)) - (6144*c4beta*pow6(M1))/(pow2(c2beta)*
        pow3(M1 - Mu)*pow3(M1 + Mu)) - (6144*s6beta*pow6(M1))/(pow2(c2beta)*
        pow3(M1 - Mu)*pow3(M1 + Mu)) - (12288*cbeta*sbeta*pow6(M1))/(pow2(
        c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) + (30720*c4beta*pow6(Mu))/(pow2(
        c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) + (6144*s6beta*pow6(Mu))/(pow2(
        c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) + (135168*pow6(Mu))/(pow2(c2beta)*
        pow3(-M1 + Mu)*pow3(M1 + Mu)) + (184320*cbeta*sbeta*pow6(Mu))/(pow2(
        c2beta)*pow3(-M1 + Mu)*pow3(M1 + Mu)) - (12*log(MR2/pow2(M1))*pow2(M1)*
        pow2(sw)*(4*(931 - 129*c4beta + 1244*c2beta*t2beta)*pow3(M1)*pow3(Mu) +
        Mu2*(277 - 743*c4beta + 1540*c2beta*t2beta)*pow4(M1) + (7409 + 1445*
        c4beta + 6324*c2beta*t2beta)*pow2(M1)*pow4(Mu) - 2*Mu*(-317 + 31*c4beta
        + 348*c2beta*t2beta)*pow5(M1) + 2*M1*(893 - 735*c4beta + 512*s6beta +
        2468*c2beta*t2beta)*pow5(Mu) + (-137 - 29*c4beta + 108*c2beta*t2beta)*
        pow6(M1) + (643 - 673*c4beta + 220*c2beta*t2beta)*pow6(Mu)))/(Mu2*pow2(
        c2beta)*pow3(M1 - Mu)*pow3(M1 + Mu)) + (840*pow2(1/cos(w))*pow7(M1))/(
        Mu2*pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) - (12*cos(4*beta - 2*w)*
        pow2(1/cos(w))*pow7(M1))/(Mu2*pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) -
        (12*cos(2*(2*beta + w))*pow2(1/cos(w))*pow7(M1))/(Mu2*pow2(c2beta)*pow2(
        M1 + Mu)*pow3(M1 - Mu)) + (864*t2beta*pow2(1/cos(w))*pow7(M1))/(c2beta*
        Mu2*pow2(M1 + Mu)*pow3(-M1 + Mu)) + (24*c4beta*pow2(1/cos(w))*pow7(M1))/(
        Mu2*pow2(c2beta)*pow2(M1 + Mu)*pow3(-M1 + Mu)) + (840*(-1 + pow2(tan(w)
        ))*pow7(M1))/(Mu2*pow2(c2beta)*pow2(M1 + Mu)*pow3(-M1 + Mu)) + (65328*
        Mu*pow2(1/cos(w))*pow4(M1)*sin(2*(beta - w)))/(pow2(c2beta)*pow2(M1 + Mu)
        *pow3(M1 - Mu)) + (3780*pow2(1/cos(w))*pow5(Mu)*sin(2*(beta - w)))/(pow2(
        c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (1536*pow2(1/cos(w))*pow6(M1)*sin(
        2*(beta - w)))/(Mu*pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (1860*
        Mu2*pow2(1/cos(w))*pow3(M1)*sin(2*(-beta + w)))/(pow2(c2beta)*pow2(M1 +
        Mu)*pow3(M1 - Mu)) + (71028*pow2(M1)*pow2(1/cos(w))*pow3(Mu)*sin(2*(-beta
        + w)))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (7308*M1*pow2(1/cos(
        w))*pow4(Mu)*sin(2*(-beta + w)))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 -
        Mu)) + (3072*pow2(1/cos(w))*pow5(M1)*sin(2*(-beta + w)))/(pow2(c2beta)*
        pow2(M1 + Mu)*pow3(M1 - Mu)) + (432*pow2(1/cos(w))*pow7(M1)*sin(2*(-beta
        + w)))/(Mu2*pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) - (1860*Mu2*pow2(
        1/cos(w))*pow3(M1)*sin(2*(beta + w)))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1
        - Mu)) - (71028*pow2(M1)*pow2(1/cos(w))*pow3(Mu)*sin(2*(beta + w)))/(
        pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (65328*Mu*pow2(1/cos(w))*
        pow4(M1)*sin(2*(beta + w)))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu))
        - (7308*M1*pow2(1/cos(w))*pow4(Mu)*sin(2*(beta + w)))/(pow2(c2beta)*pow2(
        M1 + Mu)*pow3(M1 - Mu)) - (3072*pow2(1/cos(w))*pow5(M1)*sin(2*(beta + w))
        )/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (3780*pow2(1/cos(w))*pow5(
        Mu)*sin(2*(beta + w)))/(pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (
        1536*pow2(1/cos(w))*pow6(M1)*sin(2*(beta + w)))/(Mu*pow2(c2beta)*pow2(M1
        + Mu)*pow3(M1 - Mu)) - (432*pow2(1/cos(w))*pow7(M1)*sin(2*(beta + w)))/(
        Mu2*pow2(c2beta)*pow2(M1 + Mu)*pow3(M1 - Mu)) + (pow2(cw)*pow2(1/sin(
        beta))*pow2(1/cos(w))*(179988*pow4(M1)*pow4(Mu) - 78396*c4beta*pow4(M1)*
        pow4(Mu) - 16488*s6beta*pow4(M1)*pow4(Mu) - 589824*cbeta*sbeta*pow4(M1)
        *pow4(Mu) + 460776*c2beta*t2beta*pow4(M1)*pow4(Mu) + 38358*cos(2*(2*
        beta + w))*pow4(M1)*pow4(Mu) + 49152*pow2(c2beta)*pow4(M1)*pow4(Mu) +
        47280*pow3(Mu)*pow5(M1) - 17424*c4beta*pow3(Mu)*pow5(M1) + 6144*s6beta*
        pow3(Mu)*pow5(M1) + 602112*cbeta*sbeta*pow3(Mu)*pow5(M1) - 251328*
        c2beta*t2beta*pow3(Mu)*pow5(M1) + 65544*cos(2*(2*beta + w))*pow3(Mu)*
        pow5(M1) - 111144*pow3(M1)*pow5(Mu) + 38520*c4beta*pow3(M1)*pow5(Mu) -
        864*s6beta*pow3(M1)*pow5(Mu) - 106368*c2beta*t2beta*pow3(M1)*pow5(Mu) -
        50220*cos(2*(2*beta + w))*pow3(M1)*pow5(Mu) - 46512*Mu2*pow6(M1) +
        31952*c4beta*Mu2*pow6(M1) + 5352*Mu2*s6beta*pow6(M1) + 380928*cbeta*
        Mu2*sbeta*pow6(M1) - 228248*c2beta*Mu2*t2beta*pow6(M1) + 18312*Mu2*cos(
        2*(2*beta + w))*pow6(M1) - 16384*Mu2*pow2(c2beta)*pow6(M1) - 195792*
        pow2(M1)*pow6(Mu) + 84144*c4beta*pow2(M1)*pow6(Mu) + 16056*s6beta*pow2(
        M1)*pow6(Mu) + 393216*cbeta*sbeta*pow2(M1)*pow6(Mu) - 403656*c2beta*
        t2beta*pow2(M1)*pow6(Mu) - 60936*cos(2*(2*beta + w))*pow2(M1)*pow6(Mu)
        - 49152*pow2(c2beta)*pow2(M1)*pow6(Mu) + 2928*Mu*pow7(M1) - 1488*
        c4beta*Mu*pow7(M1) - 12288*cbeta*Mu*sbeta*pow7(M1) + 1728*c2beta*Mu*
        t2beta*pow7(M1) + 360*Mu*cos(2*(2*beta + w))*pow7(M1) + 100104*M1*pow7(
        Mu) - 32664*c4beta*M1*pow7(Mu) - 5280*M1*s6beta*pow7(Mu) + 184320*
        cbeta*M1*sbeta*pow7(Mu) + 21120*c2beta*M1*t2beta*pow7(Mu) - 35268*M1*
        cos(2*(2*beta + w))*pow7(Mu) - 1680*pow8(M1) + 48*c4beta*pow8(M1) +
        1728*c2beta*t2beta*pow8(M1) - 408*cos(2*(2*beta + w))*pow8(M1) +
        103164*pow8(Mu) - 50804*c4beta*pow8(Mu) - 4920*s6beta*pow8(Mu) -
        196608*cbeta*sbeta*pow8(Mu) + 227768*c2beta*t2beta*pow8(Mu) - 14910*
        cos(2*(2*beta + w))*pow8(Mu) + 16384*pow2(c2beta)*pow8(Mu) - 6*cos(4*
        beta - 2*w)*(-6393*pow4(M1)*pow4(Mu) - 10924*pow3(Mu)*pow5(M1) + 8370*
        pow3(M1)*pow5(Mu) - 3052*Mu2*pow6(M1) + 10156*pow2(M1)*pow6(Mu) - 60*
        Mu*pow7(M1) + 5878*M1*pow7(Mu) + 68*pow8(M1) + 2485*pow8(Mu)) - 4*cos(
        2*w)*(-3*(4369 - 2112*c4beta - 49152*cbeta*sbeta + 30464*c2beta*t2beta
        + 4096*pow2(c2beta))*pow4(M1)*pow4(Mu) + 12*(1573 + 512*c4beta - 128*
        s6beta - 12544*cbeta*sbeta)*pow3(Mu)*pow5(M1) + 6*(-1247 + 16384*
        c2beta*t2beta)*pow3(M1)*pow5(Mu) + 4*Mu2*(2679 + 496*c4beta - 24576*
        cbeta*sbeta + 7616*c2beta*t2beta + 1024*pow2(c2beta))*pow6(M1) + 12*(
        851 - 1040*c4beta - 8192*cbeta*sbeta + 7616*c2beta*t2beta + 1024*pow2(
        c2beta))*pow2(M1)*pow6(Mu) - 180*Mu*pow7(M1) - 6*M1*(4589 + 1024*c4beta
        - 256*s6beta + 7680*cbeta*sbeta)*pow7(Mu) + 204*pow8(M1) + (-24345 +
        4160*c4beta + 49152*cbeta*sbeta - 30464*c2beta*t2beta - 4096*pow2(
        c2beta))*pow8(Mu)) - 70428*pow4(M1)*pow4(Mu)*sin(2*beta) + 46356*
        c4beta*pow4(M1)*pow4(Mu)*sin(2*beta) - 589824*cbeta*sbeta*pow4(M1)*
        pow4(Mu)*sin(2*beta) + 99792*c2beta*t2beta*pow4(M1)*pow4(Mu)*sin(2*
        beta) + 49152*pow2(c2beta)*pow4(M1)*pow4(Mu)*sin(2*beta) + 47280*pow3(
        Mu)*pow5(M1)*sin(2*beta) - 17424*c4beta*pow3(Mu)*pow5(M1)*sin(2*beta) +
        6144*s6beta*pow3(Mu)*pow5(M1)*sin(2*beta) + 602112*cbeta*sbeta*pow3(Mu)
        *pow5(M1)*sin(2*beta) - 251328*c2beta*t2beta*pow3(Mu)*pow5(M1)*sin(2*
        beta) - 120936*pow3(M1)*pow5(Mu)*sin(2*beta) + 44088*c4beta*pow3(M1)*
        pow5(Mu)*sin(2*beta) - 120864*c2beta*t2beta*pow3(M1)*pow5(Mu)*sin(2*
        beta) + 35328*Mu2*pow6(M1)*sin(2*beta) - 8704*c4beta*Mu2*pow6(M1)*sin(
        2*beta) + 380928*cbeta*Mu2*sbeta*pow6(M1)*sin(2*beta) - 110336*c2beta*
        Mu2*t2beta*pow6(M1)*sin(2*beta) - 16384*Mu2*pow2(c2beta)*pow6(M1)*sin(
        2*beta) + 49728*pow2(M1)*pow6(Mu)*sin(2*beta) - 37824*c4beta*pow2(M1)*
        pow6(Mu)*sin(2*beta) + 393216*cbeta*sbeta*pow2(M1)*pow6(Mu)*sin(2*beta)
        - 49920*c2beta*t2beta*pow2(M1)*pow6(Mu)*sin(2*beta) - 49152*pow2(
        c2beta)*pow2(M1)*pow6(Mu)*sin(2*beta) + 2928*Mu*pow7(M1)*sin(2*beta) -
        1488*c4beta*Mu*pow7(M1)*sin(2*beta) - 12288*cbeta*Mu*sbeta*pow7(M1)*
        sin(2*beta) + 1728*c2beta*Mu*t2beta*pow7(M1)*sin(2*beta) + 109896*M1*
        pow7(Mu)*sin(2*beta) - 38232*c4beta*M1*pow7(Mu)*sin(2*beta) - 6144*M1*
        s6beta*pow7(Mu)*sin(2*beta) + 184320*cbeta*M1*sbeta*pow7(Mu)*sin(2*
        beta) + 35616*c2beta*M1*t2beta*pow7(Mu)*sin(2*beta) - 1680*pow8(M1)*
        sin(2*beta) + 48*c4beta*pow8(M1)*sin(2*beta) + 1728*c2beta*t2beta*pow8(
        M1)*sin(2*beta) + 26220*pow8(Mu)*sin(2*beta) - 12932*c4beta*pow8(Mu)*
        sin(2*beta) - 196608*cbeta*sbeta*pow8(Mu)*sin(2*beta) + 117104*c2beta*
        t2beta*pow8(Mu)*sin(2*beta) + 16384*pow2(c2beta)*pow8(Mu)*sin(2*beta) +
        9429*pow4(M1)*pow4(Mu)*sin(6*beta - 2*w) + 1644*pow3(Mu)*pow5(M1)*sin(
        6*beta - 2*w) + 10374*pow3(M1)*pow5(Mu)*sin(6*beta - 2*w) + 2436*Mu2*
        pow6(M1)*sin(6*beta - 2*w) - 13764*pow2(M1)*pow6(Mu)*sin(6*beta - 2*w)
        - 372*Mu*pow7(M1)*sin(6*beta - 2*w) - 14910*M1*pow7(Mu)*sin(6*beta - 2*
        w) + 12*pow8(M1)*sin(6*beta - 2*w) - 1377*pow8(Mu)*sin(6*beta - 2*w) -
        1563*pow4(M1)*pow4(Mu)*sin(2*(beta - w)) - 12672*c4beta*pow4(M1)*pow4(
        Mu)*sin(2*(beta - w)) - 294912*cbeta*sbeta*pow4(M1)*pow4(Mu)*sin(2*(
        beta - w)) + 182784*c2beta*t2beta*pow4(M1)*pow4(Mu)*sin(2*(beta - w)) +
        24576*pow2(c2beta)*pow4(M1)*pow4(Mu)*sin(2*(beta - w)) - 101652*pow3(
        Mu)*pow5(M1)*sin(2*(beta - w)) - 12288*c4beta*pow3(Mu)*pow5(M1)*sin(2*(
        beta - w)) + 3072*s6beta*pow3(Mu)*pow5(M1)*sin(2*(beta - w)) + 301056*
        cbeta*sbeta*pow3(Mu)*pow5(M1)*sin(2*(beta - w)) + 75558*pow3(M1)*pow5(
        Mu)*sin(2*(beta - w)) - 196608*c2beta*t2beta*pow3(M1)*pow5(Mu)*sin(2*(
        beta - w)) - 37692*Mu2*pow6(M1)*sin(2*(beta - w)) - 3968*c4beta*Mu2*
        pow6(M1)*sin(2*(beta - w)) + 196608*cbeta*Mu2*sbeta*pow6(M1)*sin(2*(
        beta - w)) - 60928*c2beta*Mu2*t2beta*pow6(M1)*sin(2*(beta - w)) - 8192*
        Mu2*pow2(c2beta)*pow6(M1)*sin(2*(beta - w)) + 25596*pow2(M1)*pow6(Mu)*
        sin(2*(beta - w)) + 24960*c4beta*pow2(M1)*pow6(Mu)*sin(2*(beta - w)) +
        196608*cbeta*sbeta*pow2(M1)*pow6(Mu)*sin(2*(beta - w)) - 182784*c2beta*
        t2beta*pow2(M1)*pow6(Mu)*sin(2*(beta - w)) - 24576*pow2(c2beta)*pow2(
        M1)*pow6(Mu)*sin(2*(beta - w)) - 372*Mu*pow7(M1)*sin(2*(beta - w)) +
        75426*M1*pow7(Mu)*sin(2*(beta - w)) + 12288*c4beta*M1*pow7(Mu)*sin(2*(
        beta - w)) - 3072*M1*s6beta*pow7(Mu)*sin(2*(beta - w)) + 92160*cbeta*
        M1*sbeta*pow7(Mu)*sin(2*(beta - w)) + 12*pow8(M1)*sin(2*(beta - w)) +
        62607*pow8(Mu)*sin(2*(beta - w)) - 8320*c4beta*pow8(Mu)*sin(2*(beta -
        w)) - 98304*cbeta*sbeta*pow8(Mu)*sin(2*(beta - w)) + 60928*c2beta*
        t2beta*pow8(Mu)*sin(2*(beta - w)) + 8192*pow2(c2beta)*pow8(Mu)*sin(2*(
        beta - w)) - 1563*pow4(M1)*pow4(Mu)*sin(2*(beta + w)) - 12672*c4beta*
        pow4(M1)*pow4(Mu)*sin(2*(beta + w)) - 294912*cbeta*sbeta*pow4(M1)*pow4(
        Mu)*sin(2*(beta + w)) + 182784*c2beta*t2beta*pow4(M1)*pow4(Mu)*sin(2*(
        beta + w)) + 24576*pow2(c2beta)*pow4(M1)*pow4(Mu)*sin(2*(beta + w)) -
        101652*pow3(Mu)*pow5(M1)*sin(2*(beta + w)) - 12288*c4beta*pow3(Mu)*
        pow5(M1)*sin(2*(beta + w)) + 3072*s6beta*pow3(Mu)*pow5(M1)*sin(2*(beta
        + w)) + 301056*cbeta*sbeta*pow3(Mu)*pow5(M1)*sin(2*(beta + w)) + 75558*
        pow3(M1)*pow5(Mu)*sin(2*(beta + w)) - 196608*c2beta*t2beta*pow3(M1)*
        pow5(Mu)*sin(2*(beta + w)) - 37692*Mu2*pow6(M1)*sin(2*(beta + w)) -
        3968*c4beta*Mu2*pow6(M1)*sin(2*(beta + w)) + 196608*cbeta*Mu2*sbeta*
        pow6(M1)*sin(2*(beta + w)) - 60928*c2beta*Mu2*t2beta*pow6(M1)*sin(2*(
        beta + w)) - 8192*Mu2*pow2(c2beta)*pow6(M1)*sin(2*(beta + w)) + 25596*
        pow2(M1)*pow6(Mu)*sin(2*(beta + w)) + 24960*c4beta*pow2(M1)*pow6(Mu)*
        sin(2*(beta + w)) + 196608*cbeta*sbeta*pow2(M1)*pow6(Mu)*sin(2*(beta +
        w)) - 182784*c2beta*t2beta*pow2(M1)*pow6(Mu)*sin(2*(beta + w)) - 24576*
        pow2(c2beta)*pow2(M1)*pow6(Mu)*sin(2*(beta + w)) - 372*Mu*pow7(M1)*sin(
        2*(beta + w)) + 75426*M1*pow7(Mu)*sin(2*(beta + w)) + 12288*c4beta*M1*
        pow7(Mu)*sin(2*(beta + w)) - 3072*M1*s6beta*pow7(Mu)*sin(2*(beta + w))
        + 92160*cbeta*M1*sbeta*pow7(Mu)*sin(2*(beta + w)) + 12*pow8(M1)*sin(2*(
        beta + w)) + 62607*pow8(Mu)*sin(2*(beta + w)) - 8320*c4beta*pow8(Mu)*
        sin(2*(beta + w)) - 98304*cbeta*sbeta*pow8(Mu)*sin(2*(beta + w)) +
        60928*c2beta*t2beta*pow8(Mu)*sin(2*(beta + w)) + 8192*pow2(c2beta)*
        pow8(Mu)*sin(2*(beta + w)) + 9429*pow4(M1)*pow4(Mu)*sin(2*(3*beta + w))
        + 1644*pow3(Mu)*pow5(M1)*sin(2*(3*beta + w)) + 10374*pow3(M1)*pow5(Mu)*
        sin(2*(3*beta + w)) + 2436*Mu2*pow6(M1)*sin(2*(3*beta + w)) - 13764*
        pow2(M1)*pow6(Mu)*sin(2*(3*beta + w)) - 372*Mu*pow7(M1)*sin(2*(3*beta +
        w)) - 14910*M1*pow7(Mu)*sin(2*(3*beta + w)) + 12*pow8(M1)*sin(2*(3*beta
        + w)) - 1377*pow8(Mu)*sin(2*(3*beta + w))))/(2.*Mu2*pow2(c2beta)*pow2(1
        + 1/tan(beta))*pow3(M1 - Mu)*pow3(M1 + Mu)) + (3*log(MR2/Mu2)*(4*(-33116*
        pow2(cw)*pow3(M1)*pow3(Mu) + 16384*pow2(sbeta)*pow3(M1)*pow3(Mu) -
        16384*c4beta*pow2(sbeta)*pow3(M1)*pow3(Mu) + 4096*s6beta*pow2(sbeta)*
        pow3(M1)*pow3(Mu) + 4096*c2beta*t2beta*pow2(sbeta)*pow3(M1)*pow3(Mu) -
        16384*pow2(cw)*pow2(sbeta)*pow3(M1)*pow3(Mu) + 16384*c4beta*pow2(cw)*
        pow2(sbeta)*pow3(M1)*pow3(Mu) - 4096*s6beta*pow2(cw)*pow2(sbeta)*pow3(
        M1)*pow3(Mu) - 4096*c2beta*t2beta*pow2(cw)*pow2(sbeta)*pow3(M1)*pow3(
        Mu) - 111307*Mu2*pow4(M1) + 77824*Mu2*pow2(sbeta)*pow4(M1) - 53248*
        c4beta*Mu2*pow2(sbeta)*pow4(M1) + 65536*c2beta*Mu2*t2beta*pow2(sbeta)*
        pow4(M1) - 58624*Mu2*pow2(cw)*pow2(sbeta)*pow4(M1) + 14080*c4beta*Mu2*
        pow2(cw)*pow2(sbeta)*pow4(M1) - 74752*c2beta*Mu2*t2beta*pow2(cw)*pow2(
        sbeta)*pow4(M1) - 54333*pow2(cw)*pow2(M1)*pow4(Mu) - 49152*pow2(M1)*
        pow2(sbeta)*pow4(Mu) + 49152*c4beta*pow2(M1)*pow2(sbeta)*pow4(Mu) +
        29952*pow2(cw)*pow2(M1)*pow2(sbeta)*pow4(Mu) - 9984*c4beta*pow2(cw)*
        pow2(M1)*pow2(sbeta)*pow4(Mu) + 9216*c2beta*t2beta*pow2(cw)*pow2(M1)*
        pow2(sbeta)*pow4(Mu) - 49550*Mu*pow5(M1) + 32768*Mu*pow2(sbeta)*pow5(
        M1) + 65536*c2beta*Mu*t2beta*pow2(sbeta)*pow5(M1) - 32768*Mu*pow2(cw)*
        pow2(sbeta)*pow5(M1) - 65536*c2beta*Mu*t2beta*pow2(cw)*pow2(sbeta)*
        pow5(M1) - 29646*M1*pow5(Mu) + Mu*pow2(tan(w))*(33116*Mu2*pow2(cw)*
        pow3(M1) + 54333*pow2(cw)*pow2(M1)*pow3(Mu) + 111307*Mu*pow4(M1) +
        29646*M1*pow4(Mu) + 49550*pow5(M1) + 47949*pow5(Mu)) + 12288*pow2(
        sbeta)*pow6(M1) + 12288*c4beta*pow2(sbeta)*pow6(M1) - 18688*pow2(cw)*
        pow2(sbeta)*pow6(M1) + 768*c4beta*pow2(cw)*pow2(sbeta)*pow6(M1) + 3072*
        c2beta*t2beta*pow2(cw)*pow2(sbeta)*pow6(M1) - 47949*pow6(Mu) + 16384*
        pow2(sbeta)*pow6(Mu) - 16384*c4beta*pow2(sbeta)*pow6(Mu) - 9984*pow2(
        cw)*pow2(sbeta)*pow6(Mu) + 3328*c4beta*pow2(cw)*pow2(sbeta)*pow6(Mu) -
        3072*c2beta*t2beta*pow2(cw)*pow2(sbeta)*pow6(Mu) + 512*cbeta*sbeta*(16*
        (-4 + 4*c4beta - s6beta - c2beta*t2beta)*(-1 + pow2(cw))*pow3(M1)*pow3(
        Mu) + Mu2*(304 + 256*c2beta*t2beta - (229 + 292*c2beta*t2beta)*pow2(cw)
        + c4beta*(-208 + 55*pow2(cw)))*pow4(M1) - 3*(64 - 3*(13 + 4*c2beta*
        t2beta)*pow2(cw) + c4beta*(-64 + 13*pow2(cw)))*pow2(M1)*pow4(Mu) - 128*
        Mu*(1 + 2*c2beta*t2beta)*(-1 + pow2(cw))*pow5(M1) + (48 + (-73 + 12*
        c2beta*t2beta)*pow2(cw) + 3*c4beta*(16 + pow2(cw)))*pow6(M1) + (64 - 3*
        (13 + 4*c2beta*t2beta)*pow2(cw) + c4beta*(-64 + 13*pow2(cw)))*pow6(Mu))
        + 256*pow2(cbeta)*(16*(-4 + 4*c4beta - s6beta - c2beta*t2beta)*(-1 +
        pow2(cw))*pow3(M1)*pow3(Mu) + Mu2*(304 + 256*c2beta*t2beta - (229 +
        292*c2beta*t2beta)*pow2(cw) + c4beta*(-208 + 55*pow2(cw)))*pow4(M1) -
        3*(64 - 3*(13 + 4*c2beta*t2beta)*pow2(cw) + c4beta*(-64 + 13*pow2(cw)))
        *pow2(M1)*pow4(Mu) - 128*Mu*(1 + 2*c2beta*t2beta)*(-1 + pow2(cw))*pow5(
        M1) + (48 + (-73 + 12*c2beta*t2beta)*pow2(cw) + 3*c4beta*(16 + pow2(cw)
        ))*pow6(M1) + (64 - 3*(13 + 4*c2beta*t2beta)*pow2(cw) + c4beta*(-64 +
        13*pow2(cw)))*pow6(Mu))) + pow2(1/cos(w))*(-2*(M1 + Mu)*cos(4*beta - 2*w)
        *(2*Mu2*(-1713 + 421*pow2(cw))*pow3(M1) - 34*(-595 + 519*pow2(cw))*
        pow2(M1)*pow3(Mu) + Mu*(-27827 + 26535*pow2(cw))*pow4(M1) + M1*(-4831 +
        6123*pow2(cw))*pow4(Mu) + (577 + 715*pow2(cw))*pow5(M1) + (-12371 +
        11079*pow2(cw))*pow5(Mu)) + 4*cos(2*w)*(33116*pow3(M1)*pow3(Mu) +
        59599*Mu2*pow2(cw)*pow4(M1) + 106041*pow2(M1)*pow4(Mu) + 49550*Mu*pow2(
        cw)*pow5(M1) + 29646*M1*pow2(cw)*pow5(Mu) + (7135 + 10101*pow2(cw))*
        pow6(M1) + 30713*pow2(cw)*pow6(Mu)) + (M1 + Mu)*(-102376*Mu2*pow3(M1) +
        1816*c4beta*Mu2*pow3(M1) + 3220*Mu2*s6beta*pow3(M1) - 100460*c2beta*
        Mu2*t2beta*pow3(M1) + 57352*Mu2*pow2(cw)*pow3(M1) + 2568*c4beta*Mu2*
        pow2(cw)*pow3(M1) - 3332*Mu2*s6beta*pow2(cw)*pow3(M1) + 51452*c2beta*
        Mu2*t2beta*pow2(cw)*pow3(M1) + 155192*pow2(M1)*pow3(Mu) - 27848*c4beta*
        pow2(M1)*pow3(Mu) - 2812*s6beta*pow2(M1)*pow3(Mu) + 179716*c2beta*
        t2beta*pow2(M1)*pow3(Mu) - 110168*pow2(cw)*pow2(M1)*pow3(Mu) + 23464*
        c4beta*pow2(cw)*pow2(M1)*pow3(Mu) + 2924*s6beta*pow2(cw)*pow2(M1)*pow3(
        Mu) - 130708*c2beta*t2beta*pow2(cw)*pow2(M1)*pow3(Mu) - 91420*Mu*pow4(
        M1) + 37732*c4beta*Mu*pow4(M1) + 798*Mu*s6beta*pow4(M1) - 128098*
        c2beta*Mu*t2beta*pow4(M1) + 68908*Mu*pow2(cw)*pow4(M1) - 35540*c4beta*
        Mu*pow2(cw)*pow4(M1) - 854*Mu*s6beta*pow2(cw)*pow4(M1) + 103594*c2beta*
        Mu*t2beta*pow2(cw)*pow4(M1) + 38644*M1*pow4(Mu) + 10356*c4beta*M1*pow4(
        Mu) + 3478*M1*s6beta*pow4(Mu) + 31510*c2beta*M1*t2beta*pow4(Mu) -
        16132*M1*pow2(cw)*pow4(Mu) - 12548*c4beta*M1*pow2(cw)*pow4(Mu) - 3422*
        M1*s6beta*pow2(cw)*pow4(Mu) - 7006*c2beta*M1*t2beta*pow2(cw)*pow4(Mu) +
        24820*pow5(M1) - 1932*c4beta*pow5(M1) + 470*s6beta*pow5(M1) + 26966*
        c2beta*t2beta*pow5(M1) - 2308*pow2(cw)*pow5(M1) - 260*c4beta*pow2(cw)*
        pow5(M1) - 414*s6beta*pow2(cw)*pow5(M1) - 2462*c2beta*t2beta*pow2(cw)*
        pow5(M1) - 86300*pow5(Mu) + 16740*c4beta*pow5(Mu) + 990*s6beta*pow5(Mu)
        - 101794*c2beta*t2beta*pow5(Mu) + 63788*pow2(cw)*pow5(Mu) - 14548*
        c4beta*pow2(cw)*pow5(Mu) - 1046*s6beta*pow2(cw)*pow5(Mu) + 77290*
        c2beta*t2beta*pow2(cw)*pow5(Mu) - 2*cos(2*(2*beta + w))*(2*Mu2*(-1713 +
        421*pow2(cw))*pow3(M1) - 34*(-595 + 519*pow2(cw))*pow2(M1)*pow3(Mu) +
        Mu*(-27827 + 26535*pow2(cw))*pow4(M1) + M1*(-4831 + 6123*pow2(cw))*
        pow4(Mu) + (577 + 715*pow2(cw))*pow5(M1) + (-12371 + 11079*pow2(cw))*
        pow5(Mu)) + (2*Mu2*(4087 - 4587*pow2(cw))*pow3(M1) + 2*(-3701 + 4201*
        pow2(cw))*pow2(M1)*pow3(Mu) + Mu*(2069 - 2569*pow2(cw))*pow4(M1) + M1*(
        2217 - 1717*pow2(cw))*pow4(Mu) + (361 + 139*pow2(cw))*pow5(M1) + (3797
        - 4297*pow2(cw))*pow5(Mu))*sin(6*beta - 2*w) - 107666*Mu2*pow3(M1)*sin(
        2*(beta - w)) + 33322*Mu2*pow2(cw)*pow3(M1)*sin(2*(beta - w)) + 208278*
        pow2(M1)*pow3(Mu)*sin(2*(beta - w)) - 133934*pow2(cw)*pow2(M1)*pow3(Mu)
        *sin(2*(beta - w)) - 166827*Mu*pow4(M1)*sin(2*(beta - w)) + 129655*Mu*
        pow2(cw)*pow4(M1)*sin(2*(beta - w)) + 29033*M1*pow4(Mu)*sin(2*(beta -
        w)) + 8139*M1*pow2(cw)*pow4(Mu)*sin(2*(beta - w)) + 15657*pow5(M1)*sin(
        2*(beta - w)) + 21515*pow2(cw)*pow5(M1)*sin(2*(beta - w)) - 116715*
        pow5(Mu)*sin(2*(beta - w)) + 79543*pow2(cw)*pow5(Mu)*sin(2*(beta - w))
        - 107666*Mu2*pow3(M1)*sin(2*(beta + w)) + 33322*Mu2*pow2(cw)*pow3(M1)*
        sin(2*(beta + w)) + 208278*pow2(M1)*pow3(Mu)*sin(2*(beta + w)) -
        133934*pow2(cw)*pow2(M1)*pow3(Mu)*sin(2*(beta + w)) - 166827*Mu*pow4(
        M1)*sin(2*(beta + w)) + 129655*Mu*pow2(cw)*pow4(M1)*sin(2*(beta + w)) +
        29033*M1*pow4(Mu)*sin(2*(beta + w)) + 8139*M1*pow2(cw)*pow4(Mu)*sin(2*(
        beta + w)) + 15657*pow5(M1)*sin(2*(beta + w)) + 21515*pow2(cw)*pow5(M1)
        *sin(2*(beta + w)) - 116715*pow5(Mu)*sin(2*(beta + w)) + 79543*pow2(cw)
        *pow5(Mu)*sin(2*(beta + w)) + 8174*Mu2*pow3(M1)*sin(2*(3*beta + w)) -
        9174*Mu2*pow2(cw)*pow3(M1)*sin(2*(3*beta + w)) - 7402*pow2(M1)*pow3(Mu)
        *sin(2*(3*beta + w)) + 8402*pow2(cw)*pow2(M1)*pow3(Mu)*sin(2*(3*beta +
        w)) + 2069*Mu*pow4(M1)*sin(2*(3*beta + w)) - 2569*Mu*pow2(cw)*pow4(M1)*
        sin(2*(3*beta + w)) + 2217*M1*pow4(Mu)*sin(2*(3*beta + w)) - 1717*M1*
        pow2(cw)*pow4(Mu)*sin(2*(3*beta + w)) + 361*pow5(M1)*sin(2*(3*beta + w)
        ) + 139*pow2(cw)*pow5(M1)*sin(2*(3*beta + w)) + 3797*pow5(Mu)*sin(2*(3*
        beta + w)) - 4297*pow2(cw)*pow5(Mu)*sin(2*(3*beta + w))))))/(4.*pow2(
        c2beta)*pow2(cbeta + sbeta)*pow3(M1 - Mu)*pow3(M1 + Mu)))/98304.;

      return lim;
   }

   const double exact =
      (2*lmL3MR*pow2(cw) + 6*lmQ3MR*pow2(cw) + 8*lmUMR*pow2(cw) + 2*log(pow2(p.MA)/
        MR2)*pow2(cw) + 16*log(pow2(M2)/MR2)*pow2(cw) + 2*log(pow2(mL1)/MR2)*
        pow2(cw) + 2*log(pow2(mL2)/MR2)*pow2(cw) + 6*log(pow2(mQ1)/MR2)*pow2(
        cw) + 6*log(pow2(mQ2)/MR2)*pow2(cw) + (3*log(MR2/pow2(M1))*pow2(sw)*
        pow3(M1)*(-((-Mu2 + pow2(M2))*(-((1 + c4beta)*M2*Mu) - 2*Mu2*(s2beta +
        s6beta) + 2*(s2beta + s6beta)*pow2(M2))*pow3(Mu)) - 2*Mu2*pow2(M1)*((5
        - 3*c4beta)*M2*Mu2 + 4*Mu*s2beta*pow2(M2) + (1 + c4beta)*pow3(M2) + 4*
        s2beta*pow3(Mu)) + ((11 - 5*c4beta)*M2*Mu2 + 8*Mu*s2beta*pow2(M2) + (1
        + c4beta)*pow3(M2) + 8*s2beta*pow3(Mu))*pow4(M1) + 2*(1 + c4beta)*pow3(
        M1)*(-3*Mu2*pow2(M2) + pow4(M2) + 2*pow4(Mu)) - 2*M1*Mu2*(-((5 + 7*
        c4beta)*Mu2*pow2(M2)) + 4*M2*s2beta*pow3(Mu) + 3*(1 + c4beta)*pow4(M2)
        + (5 + 3*c4beta)*pow4(Mu)) + 2*(-((-1 + c4beta)*Mu2) + 4*M2*Mu*s2beta +
        2*pow2(M2))*pow5(M1)))/(pow2(c2beta)*pow2(M2 - Mu)*pow2(M2 + Mu)*pow3(
        M1 - Mu)*pow3(M1 + Mu)) + (3*log(MR2/pow2(M2))*pow3(M2)*((1 + c4beta)*(
        -1 + pow2(cw))*pow2(-Mu2 + pow2(M2))*pow3(M1) + M1*Mu*(-1 + pow2(cw))*(
        -Mu2 + pow2(M2))*(8*M2*Mu2*s2beta + (11 - 5*c4beta)*Mu*pow2(M2) + 8*
        s2beta*pow3(M2) + (1 + c4beta)*pow3(Mu)) + (M2*Mu2*(18 + c4beta*(18 -
        13*pow2(cw)) + 35*pow2(cw)) + 48*Mu*s2beta*pow2(cw)*pow2(M2) - (6 - 13*
        pow2(cw) + 3*c4beta*(2 + pow2(cw)))*pow3(M2) + 2*(3*s6beta + s2beta*(3
        + 8*pow2(cw)))*pow3(Mu))*pow4(M1) + 2*pow2(M1)*(Mu2*(7 - 14*pow2(cw) +
        c4beta*(7 + 2*pow2(cw)))*pow3(M2) + 4*s2beta*(1 - 13*pow2(cw))*pow2(M2)
        *pow3(Mu) + 4*Mu*s2beta*(-1 + pow2(cw))*pow4(M2) + M2*(-17 - 36*pow2(
        cw) + c4beta*(-19 + 14*pow2(cw)))*pow4(Mu) + 2*(-1 + pow2(cw))*pow5(M2)
        - 2*(3*s6beta + s2beta*(3 + 8*pow2(cw)))*pow5(Mu)) + Mu2*(-(Mu2*(8 -
        15*pow2(cw) + c4beta*(8 + pow2(cw)))*pow3(M2)) + 8*s2beta*(1 + 5*pow2(
        cw))*pow2(M2)*pow3(Mu) + 8*Mu*s2beta*(-1 + pow2(cw))*pow4(M2) + M2*(22
        + c4beta*(18 - 13*pow2(cw)) + 31*pow2(cw))*pow4(Mu) - 2*(-1 + c4beta)*(
        -1 + pow2(cw))*pow5(M2) + 2*(3*s6beta + s2beta*(3 + 8*pow2(cw)))*pow5(
        Mu))))/(pow2(c2beta)*pow2(M1 - Mu)*pow2(M1 + Mu)*pow3(-M2 + Mu)*pow3(M2
        + Mu)) + (pow4(M1)*(2*Mu2*(-27 - 68*pow2(cw) + c4beta*(-33 + 34*pow2(
        cw)))*pow2(M2) - 9*Mu*(s2beta + s6beta + 16*s2beta*pow2(cw))*pow3(M2) -
        3*M2*(3*s6beta + s2beta*(-5 + 56*pow2(cw)))*pow3(Mu) + 2*(6 - 23*pow2(
        cw) + c4beta*(6 + pow2(cw)))*pow4(M2) + (24 - 55*pow2(cw) + c4beta*(12
        + 5*pow2(cw)))*pow4(Mu)) - 3*(-1 + pow2(cw))*((11 - 5*c4beta)*M2*Mu2 +
        8*Mu*s2beta*pow2(M2) + (1 + c4beta)*pow3(M2) + 8*s2beta*pow3(Mu))*pow5(
        M1) - 3*(-1 + pow2(cw))*pow3(M1)*(-2*(1 + c4beta)*Mu2*pow3(M2) + 2*(-
        s2beta + s6beta)*pow2(M2)*pow3(Mu) - Mu*(s2beta + s6beta)*pow4(M2) + (-
        5 + 3*c4beta)*M2*pow4(Mu) + (1 + c4beta)*pow5(M2) - (5*s2beta + s6beta)
        *pow5(Mu)) - 6*(-1 + pow2(cw))*(-((-1 + c4beta)*Mu2) + 4*M2*Mu*s2beta +
        2*pow2(M2))*pow6(M1) - 3*M1*Mu*(-1 + pow2(cw))*((-5 + 3*c4beta)*pow3(
        M2)*pow3(Mu) + Mu2*(7*s2beta - s6beta)*pow4(M2) + 2*(-9*s2beta +
        s6beta)*pow2(M2)*pow4(Mu) + (11 - 5*c4beta)*Mu*pow5(M2) + 4*(-3 +
        c4beta)*M2*pow5(Mu) + 8*s2beta*pow6(M2) - (5*s2beta + s6beta)*pow6(Mu))
        - 2*pow2(M1)*(-3*(s2beta + 3*s6beta + 50*s2beta*pow2(cw))*pow3(M2)*
        pow3(Mu) + Mu2*(15 - 7*c4beta*(-3 + pow2(cw)) - 49*pow2(cw))*pow4(M2) +
        2*(-27 - 68*pow2(cw) + c4beta*(-39 + 40*pow2(cw)))*pow2(M2)*pow4(Mu) +
        12*Mu*s2beta*(-1 + pow2(cw))*pow5(M2) - 3*M2*(3*s6beta + s2beta*(-7 +
        58*pow2(cw)))*pow5(Mu) + 6*(-1 + pow2(cw))*pow6(M2) + (27 - 58*pow2(cw)
        + c4beta*(15 + 2*pow2(cw)))*pow6(Mu)) + Mu2*(-3*(3*s6beta + s2beta*(7 +
        44*pow2(cw)))*pow3(M2)*pow3(Mu) + 2*Mu2*(12 - 29*pow2(cw) + c4beta*(6 +
        pow2(cw)))*pow4(M2) + 2*(-39 - 56*pow2(cw) + c4beta*(-27 + 28*pow2(cw))
        )*pow2(M2)*pow4(Mu) - 24*Mu*s2beta*(-1 + pow2(cw))*pow5(M2) - 3*M2*(3*
        s6beta + s2beta*(7 + 44*pow2(cw)))*pow5(Mu) + 6*(-1 + c4beta)*(-1 +
        pow2(cw))*pow6(M2) + (12 - 43*pow2(cw) + c4beta*(12 + 5*pow2(cw)))*
        pow6(Mu)))/(pow2(c2beta)*pow2(M1 - Mu)*pow2(M2 - Mu)*pow2(M1 + Mu)*
        pow2(M2 + Mu)) + (3*log(MR2/Mu2)*((Mu2*(-1 + s2beta)*s2beta*(M2*(7*M2 -
        3*Mu)*Mu2 + pow2(M1)*(5*M2*Mu + 7*Mu2 - 8*pow2(M2)) + M1*Mu*(6*M2*Mu -
        3*Mu2 + 5*pow2(M2))))/((M1 - Mu)*(-M2 + Mu)*pow2(M1 + Mu)*pow2(M2 + Mu)
        ) - (Mu2*(-1 + s2beta)*(pow2(M1)*(M2*Mu + 11*Mu2 - 16*pow2(M2)) + M1*
        Mu*(-2*M2*Mu - 7*Mu2 + pow2(M2)) + Mu2*(-7*M2*Mu - 8*Mu2 + 11*pow2(M2))
        ))/((M1 - Mu)*(-M2 + Mu)*pow2(M1 + Mu)*pow2(M2 + Mu)) - 16*pow2(s2beta)
        - (Mu2*s2beta*(M2*(7*M2 + 3*Mu)*Mu2 + pow2(M1)*(-5*M2*Mu + 7*Mu2 - 8*
        pow2(M2)) + M1*Mu*(6*M2*Mu + 3*Mu2 - 5*pow2(M2)))*pow2(cbeta + sbeta))/
        ((M1 + Mu)*(M2 + Mu)*pow2(M1 - Mu)*pow2(M2 - Mu)) + (Mu2*(Mu2*(-7*M2*Mu
        + 8*Mu2 - 11*pow2(M2)) + M1*Mu*(2*M2*Mu - 7*Mu2 + pow2(M2)) + pow2(M1)*
        (M2*Mu - 11*Mu2 + 16*pow2(M2)))*pow2(cbeta + sbeta))/((M1 + Mu)*(M2 +
        Mu)*pow2(M1 - Mu)*pow2(M2 - Mu)) - (2*pow3(M1)*(M1*Mu2*(3*(1 + c4beta)*
        M2 + Mu*(15*s2beta - s6beta)) + (1 + c4beta)*Mu2*pow2(M1) - (M2 +
        c4beta*M2 + 16*Mu*s2beta)*pow3(M1) + (-4*(-1 + c4beta)*Mu + M2*(s2beta
        + s6beta))*pow3(Mu) + (-7 + c4beta)*pow4(M1)))/((M1 - M2)*pow3(-Mu2 +
        pow2(M1))) + (2*pow3(M2)*(-5*(1 + c4beta)*Mu2*pow2(M2) - 16*Mu*s2beta*
        pow3(M2) + M2*(13*s2beta - 3*s6beta)*pow3(Mu) + M1*(9*(1 + c4beta)*M2*
        Mu2 - 3*(1 + c4beta)*pow3(M2) + 3*(s2beta + s6beta)*pow3(Mu)) + (-5 +
        3*c4beta)*pow4(M2) - 4*(-1 + c4beta)*pow4(Mu)))/((M1 - M2)*pow3(-Mu2 +
        pow2(M2))) - (2*pow2(cw)*pow3(Mu)*(M1*(-Mu2 + pow2(M2))*pow3(Mu)*(2*(-3
        + c4beta)*M2*Mu2 - 10*Mu*s2beta*pow2(M2) + (-3 + c4beta)*pow3(M2) - 2*
        s2beta*pow3(Mu)) - (-Mu2 + pow2(M2))*pow3(M1)*(-2*Mu2*(5*s2beta +
        s6beta)*pow2(M2) + (-3 + c4beta)*Mu*pow3(M2) + (-3 + c4beta)*M2*pow3(
        Mu) + (s2beta + s6beta)*pow4(M2) + (s2beta + s6beta)*pow4(Mu)) + (-Mu2
        + pow2(M2))*(-((-3 + c4beta)*M2*Mu) + 2*Mu2*s2beta + 2*s2beta*pow2(M2))
        *pow5(M1) + pow4(M1)*(-16*Mu2*s2beta*pow3(M2) + 19*(-3 + c4beta)*pow2(
        M2)*pow3(Mu) - 2*(-3 + c4beta)*Mu*pow4(M2) - 82*M2*s2beta*pow4(Mu) + 2*
        s2beta*pow5(M2) + 7*(-3 + c4beta)*pow5(Mu)) + 2*(12*M2*Mu2*s2beta - 3*(
        -3 + c4beta)*Mu*pow2(M2) + 4*s2beta*pow3(M2) - (-3 + c4beta)*pow3(Mu))*
        pow6(M1) + pow3(Mu)*(-8*s2beta*pow3(M2)*pow3(Mu) - 2*(3 + c4beta)*Mu2*
        pow4(M2) + 4*(-3 + 2*c4beta)*pow2(M2)*pow4(Mu) - 2*Mu*s2beta*pow5(M2) -
        22*M2*s2beta*pow5(Mu) + (1 + c4beta)*pow6(M2) + (-7 + c4beta)*pow6(Mu))
        + pow2(M1)*(2*(3 + 5*c4beta)*pow3(Mu)*pow4(M2) + 16*s2beta*pow3(M2)*
        pow4(Mu) - 9*(-5 + 3*c4beta)*pow2(M2)*pow5(Mu) - 3*(1 + c4beta)*Mu*
        pow6(M2) + 80*M2*s2beta*pow6(Mu) - 4*(-6 + c4beta)*pow7(Mu))))/(pow3(M1
        - Mu)*pow3(M1 + Mu)*pow3(-M2 + Mu)*pow3(M2 + Mu))))/pow2(c2beta))/48.;

   return exact;
}

double ThresholdCalculator::getDeltaYtauYb() const
{
   const double Xt2 = pow2(p.Au(2,2) - p.mu*p.vd/p.vu);
   const double mQ3 = std::sqrt(p.mq2(2,2));
   const double mU3 = std::sqrt(p.mu2(2,2));
   const double beta = std::atan(p.vu/p.vd);
   const double sbeta = std::sin(beta);

   return -(pow2(sbeta)*Xt2*F5(mQ3/mU3))/(4*mQ3*mU3);
}

double ThresholdCalculator::getDeltaLambdaYb4G32(int omitLogs) const
{
   using std::log;
   using himalaya::dilog;

   const double MR2 = pow2(p.scale);
   const double lMR = omitLogs*log(MR2);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double m32 = pow2(p.MG);
   const double mQ3 = std::sqrt(mQ32);
   const double mD3 = std::sqrt(mD32);
   const double m3 = std::sqrt(m32);
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double Xb2 = pow2(Xb);
   const double eps = 0.01*mQ3;
   const double lmQ3MR = log(mQ32) - lMR;
   const double lmD3MR = log(mD32) - lMR;

   if (std::abs(mQ3 - m3) < eps && std::abs(mD3 - m3) < eps && std::abs(mQ3 - mD3) < eps) {
      return (-8*(Xb*(-12*mQ32*Xb + 4*mQ3*pow2(Xb) - 24*pow3(mQ3) + pow3(Xb)) + 8*mQ3*
        lMR*(-3*mQ32*Xb - 3*mQ3*pow2(Xb) + 3*pow3(mQ3) + pow3(Xb)) -
        16*mQ3*log(mQ3)*(-3*mQ32*Xb - 3*mQ3*pow2(Xb) + 3*pow3(mQ3) + 3*lMR
        *pow3(mQ3) + pow3(Xb)) + 48*pow2(log(mQ3))*pow4(mQ3) + 12*
        pow2(lMR)*pow4(mQ3)))/(3.*pow4(mQ3));
   } else if(std::abs(mD3 - mQ3) < eps){
     return (-4*(-48*pow3(m3)*pow3(Xb)*pow4(mQ3) + 156*lmQ3MR*pow3(m3)*pow3(Xb)*pow4(
        mQ3) - 124*log(m32/MR2)*pow3(m3)*pow3(Xb)*pow4(mQ3) - 192*Xb2*pow4(m3)*
        pow4(mQ3) + 96*lmQ3MR*Xb2*pow4(m3)*pow4(mQ3) + 96*Xb2*log(m32/MR2)*
        pow4(m3)*pow4(mQ3) + 24*mQ32*pow4(m3)*pow4(Xb) + 54*lmQ3MR*mQ32*pow4(
        m3)*pow4(Xb) - 78*mQ32*log(m32/MR2)*pow4(m3)*pow4(Xb) - 20*m32*pow4(
        mQ3)*pow4(Xb) - 15*lmQ3MR*m32*pow4(mQ3)*pow4(Xb) + 39*m32*log(m32/MR2)*
        pow4(mQ3)*pow4(Xb) + 16*mQ32*pow3(Xb)*pow5(m3) - 264*lmQ3MR*mQ32*pow3(
        Xb)*pow5(m3) + 248*mQ32*log(m32/MR2)*pow3(Xb)*pow5(m3) - 96*Xb*pow4(
        mQ3)*pow5(m3) + 192*lmQ3MR*Xb*pow4(mQ3)*pow5(m3) - 96*Xb*log(m32/MR2)*
        pow4(mQ3)*pow5(m3) + m32*Xb*log(m32/mQ32)*(66*m32*mQ32*pow3(Xb) - 31*
        pow3(Xb)*pow4(m3) + 156*m3*Xb2*pow4(mQ3) - 39*pow3(Xb)*pow4(mQ3) + 24*
        pow3(m3)*(-11*mQ32*Xb2 + 8*pow4(mQ3)) + 124*Xb2*pow5(m3)) + 48*mQ32*
        Xb2*pow6(m3) - 48*mQ32*Xb2*log(m32/MR2)*pow6(m3) - 48*pow4(mQ3)*pow6(
        m3) + 48*log(m32/MR2)*pow4(mQ3)*pow6(m3) - 8*pow4(Xb)*pow6(m3) - 31*
        lmQ3MR*pow4(Xb)*pow6(m3) + 39*log(m32/MR2)*pow4(Xb)*pow6(m3) + 192*m32*
        Xb2*pow6(mQ3) - 240*lmQ3MR*m32*Xb2*pow6(mQ3) + 192*Xb*pow3(m3)*pow6(
        mQ3) - 192*lmQ3MR*Xb*pow3(m3)*pow6(mQ3) + 96*(m3 + 2*Xb)*dilog(1 - m32/
        mQ32)*pow3(m3)*pow6(mQ3) + 96*Xb*log(m32/MR2)*pow3(m3)*pow6(mQ3) - 96*
        lmQ3MR*Xb*log(m32/MR2)*pow3(m3)*pow6(mQ3) + 96*Xb*pow2(lmQ3MR)*pow3(m3)
        *pow6(mQ3) + 32*m3*pow3(Xb)*pow6(mQ3) - 16*lmQ3MR*m3*pow3(Xb)*pow6(mQ3)
        + 156*pow4(m3)*pow6(mQ3) - 144*lmQ3MR*pow4(m3)*pow6(mQ3) + 24*log(m32/
        MR2)*pow4(m3)*pow6(mQ3) - 48*lmQ3MR*log(m32/MR2)*pow4(m3)*pow6(mQ3) +
        72*pow2(lmQ3MR)*pow4(m3)*pow6(mQ3) + 4*pow4(Xb)*pow6(mQ3) - 8*lmQ3MR*
        pow4(Xb)*pow6(mQ3) + 124*lmQ3MR*pow3(Xb)*pow7(m3) - 124*log(m32/MR2)*
        pow3(Xb)*pow7(m3) - 144*m32*pow8(mQ3) + 144*lmQ3MR*m32*pow8(mQ3) - 96*
        m3*Xb*pow8(mQ3) - 48*Xb2*pow8(mQ3) + 96*lmQ3MR*Xb2*pow8(mQ3) - 48*m32*
        pow2(lmQ3MR)*pow8(mQ3) + 36*power10(mQ3) - 72*lmQ3MR*power10(mQ3) + 24*
        pow2(lmQ3MR)*power10(mQ3)))/(3.*pow2(m32 - mQ32)*pow6(mQ3));
   } else if(std::abs(m3 - mQ3) < eps){
     return 12 - 8*lmQ3MR + (8*(-13*mD32 + 13*mQ32 + 2*lmQ3MR*(3*mD32 - mQ3*(3*mQ3 +
        8*Xb)) + 8*Xb2))/(mD32 - mQ32) - (2*pow2(lmQ3MR)*(-2*mD32*mQ3*(27*mQ3 +
        40*Xb) + (43*mQ3 + 16*Xb)*pow3(mQ3) + 27*pow4(mD3)))/pow2(mD32 - mQ32)
        + (lmQ3MR*(-16*mQ32*Xb2 - 6*mD32*(43*mQ32 + 8*Xb2) - 256*mQ3*pow3(Xb) +
        129*pow4(mD3) + 177*pow4(mQ3) - 128*pow4(Xb)))/pow2(mD32 - mQ32) + (
        lmQ3MR*(-32*pow2(-(mQ3*Xb2) + pow3(mQ3)) + 2*(89*mQ32 + 56*Xb2)*pow4(
        mD3) + 2*lmQ3MR*mD32*(-2*mD32*mQ3*(11*mQ3 + 40*Xb) + (27*mQ3 + 16*Xb)*
        pow3(mQ3) + 11*pow4(mD3)) - mD32*(112*mQ32*Xb2 + 105*pow4(mQ3) + 32*
        pow4(Xb)) - 89*pow6(mD3)))/pow2(-(mD3*mQ32) + pow3(mD3)) + (4*(8*pow2(-
        (mQ3*Xb2) + pow3(mQ3)) - 2*(11*mQ32 + 32*mQ3*Xb + 8*Xb2)*pow4(mD3) - 8*
        mD32*dilog(1 - mQ32/mD32)*(-8*Xb*pow3(mQ3) - 4*mQ3*pow3(Xb) + 2*pow4(
        mQ3) + pow4(Xb)) + mD32*(32*mQ32*Xb2 + 64*Xb*pow3(mQ3) + 128*mQ3*pow3(
        Xb) - 3*pow4(mQ3) + 56*pow4(Xb)) + 17*pow6(mD3)))/pow2(-(mD3*mQ32) +
        pow3(mD3)) - (16*log(mD32/mQ32)*((9*mQ32 - 8*mQ3*Xb - 8*Xb2)*pow4(mD3)
        + 2*mQ32*(-6*mQ32*Xb2 - 4*Xb*pow3(mQ3) + 4*mQ3*pow3(Xb) + 3*pow4(mQ3) +
        5*pow4(Xb)) + 2*mD32*(10*mQ32*Xb2 + 8*Xb*pow3(mQ3) + 12*mQ3*pow3(Xb) -
        6*pow4(mQ3) + 7*pow4(Xb)) - 3*pow6(mD3) + 2*lmQ3MR*((-3*mQ32 + 4*Xb2)*
        pow4(mD3) + mD32*(-8*mQ32*Xb2 - 4*Xb*pow3(mQ3) - 4*mQ3*pow3(Xb) + 4*
        pow4(mQ3) - 2*pow4(Xb)) - 2*mQ32*(-2*mQ32*Xb2 - 2*Xb*pow3(mQ3) + 2*mQ3*
        pow3(Xb) + pow4(mQ3) + 2*pow4(Xb)) + pow6(mD3))))/pow3(mD32 - mQ32) - (
        16*pow2(log(mD32/mQ32))*(2*mD32*mQ32*(5*mQ32*Xb2 + 8*Xb*pow3(mQ3) + 4*
        mQ3*pow3(Xb) - 4*pow4(mQ3) - 3*pow4(Xb)) + 2*pow4(mD3)*(-7*mQ32*Xb2 -
        4*Xb*pow3(mQ3) - 2*mQ3*pow3(Xb) + 4*pow4(mQ3) - 2*pow4(Xb)) + pow4(mQ3)
        *(-2*mQ32*Xb2 - 8*Xb*pow3(mQ3) - 4*mQ3*pow3(Xb) + 3*pow4(mQ3) + 2*pow4(
        Xb)) + (-4*mQ32 + 6*Xb2)*pow6(mD3) + pow8(mD3)))/pow4(mD32 - mQ32);
   } else if(std::abs(m3 - mD3) < eps){
     return 12 - 8*lmQ3MR + (8*(-13*mD32 + 13*mQ32 + 2*lmQ3MR*(3*mD32 - 3*mQ32 + 8*
        mD3*Xb) - 8*Xb2))/(mD32 - mQ32) - (2*pow2(lmQ3MR)*(-54*mD32*mQ32 - 80*
        mD3*mQ32*Xb + 16*Xb*pow3(mD3) + 43*pow4(mD3) + 27*pow4(mQ3)))/pow2(mD32
        - mQ32) + (lmQ3MR*(-48*mQ32*Xb2 - 2*mD32*(129*mQ32 + 8*Xb2) - 256*mD3*
        pow3(Xb) + 177*pow4(mD3) + 129*pow4(mQ3) - 128*pow4(Xb)))/pow2(mD32 -
        mQ32) + 16*pow2(log(mD32/mQ32))*(-2.375 - (10*mD3*Xb)/(mD32 - mQ32) - (
        2*(3*mD32 - mQ32)*Xb2)/pow2(mD32 - mQ32) - (4*mD3*pow3(Xb))/pow2(mD32 -
        mQ32) + (2*mD32*(3*mD32 + mQ32)*pow4(Xb))/pow4(mD32 - mQ32)) + (log(
        mD32/mQ32)*(-352*mD32*mQ32*Xb2 - 512*mQ32*Xb*pow3(mD3) - 128*mD3*mQ32*
        pow3(Xb) - 384*pow3(mD3)*pow3(Xb) - 363*mQ32*pow4(mD3) + 176*Xb2*pow4(
        mD3) + 363*mD32*pow4(mQ3) + 256*mD3*Xb*pow4(mQ3) + 176*Xb2*pow4(mQ3) -
        288*mD32*pow4(Xb) - 96*mQ32*pow4(Xb) + 256*Xb*pow5(mD3) + 121*pow6(mD3)
        - 4*lmQ3MR*(-32*pow3(mD3)*(5*mQ32*Xb - pow3(Xb)) + (-57*mQ32 + 32*Xb2)*
        pow4(mD3) + 32*Xb2*pow4(mQ3) + 16*mD3*(-2*mQ32*pow3(Xb) + 5*Xb*pow4(
        mQ3)) + mD32*(-64*mQ32*Xb2 + 57*pow4(mQ3) - 32*pow4(Xb)) - 16*mQ32*
        pow4(Xb) + 80*Xb*pow5(mD3) + 19*pow6(mD3) - 19*pow6(mQ3)) + 2*lmD3MR*(-
        32*pow3(mD3)*(5*mQ32*Xb - 4*pow3(Xb)) - 33*mQ32*pow4(mD3) + 33*mD32*
        pow4(mQ3) + 80*mD3*Xb*pow4(mQ3) + 80*Xb*pow5(mD3) + 11*pow6(mD3) - 11*
        pow6(mQ3)) - 121*pow6(mQ3)))/pow3(mD32 - mQ32) + (4*(-64*mD3*mQ32*Xb*(
        mQ32 - 2*Xb2) + 64*mQ32*Xb*pow3(mD3) - (3*mQ32 + 16*Xb2)*pow4(mD3) -
        16*Xb2*pow4(mQ3) + 56*mQ32*pow4(Xb) - 8*mQ32*dilog(1 - mD32/mQ32)*(-8*
        Xb*pow3(mD3) - 4*mD3*pow3(Xb) + 2*pow4(mD3) + pow4(Xb)) + mD32*(32*
        mQ32*Xb2 - 22*pow4(mQ3) + 8*pow4(Xb)) + 8*pow6(mD3) + 17*pow6(mQ3)))/
        pow2(-(mD32*mQ3) + pow3(mQ3)) - (lmD3MR*((105*mQ32 - 64*Xb2)*pow4(mD3)
        - 112*Xb2*pow4(mQ3) - 2*lmQ3MR*mQ32*(-22*mD32*mQ32 - 80*mD3*mQ32*Xb +
        16*Xb*pow3(mD3) + 27*pow4(mD3) + 11*pow4(mQ3)) + 32*mQ32*pow4(Xb) +
        mD32*(112*mQ32*Xb2 - 178*pow4(mQ3) + 32*pow4(Xb)) + 32*pow6(mD3) + 89*
        pow6(mQ3)))/pow2(-(mD32*mQ3) + pow3(mQ3));
   }

   return 16*(-2*pow2(lMR) + (16*Xb*(m3*(m32 - mD32)*(m32 - mQ32)*log(mD3) - (m32 -
        mQ32)*dilog(1 - m32/mD32)*pow3(m3) + (m32 - mD32)*dilog(1 - m32/mQ32)*
        pow3(m3) + lMR*(m32 - mQ32)*log(mD3)*pow3(m3) + log(m3)*(lMR*(-mD32 +
        mQ32) + 2*(m32 - mQ32)*log(mD3) - 2*(m32 - mD32)*log(mQ3))*pow3(m3) -
        2*(m32 - mQ32)*pow2(log(mD3))*pow3(m3) + 2*(m32 - mD32)*pow2(log(mQ3))*
        pow3(m3) - (m32 - mD32)*log(mQ3)*(-(m3*mQ32) + (1 + lMR)*pow3(m3))))/((
        m32 - mD32)*(m32 - mQ32)*(mD32 - mQ32)) - (8*m3*((2*m32 - mD32 - mQ32)*
        dilog(1 - m32/mD32) + (-2*m32 + mD32 + mQ32)*dilog(1 - m32/mQ32) + 2*(-
        ((2 + lMR)*(mD32 - mQ32)) + ((3 + lMR)*mD32 + (1 + lMR)*mQ32 - 4*m32*
        log(m3))*log(mD3) - ((1 + lMR)*mD32 + (3 + lMR)*mQ32 - 4*m32*log(m3))*
        log(mQ3) + (2*m32 - mD32 - mQ32)*pow2(log(mD3)) + (-2*m32 + mD32 +
        mQ32)*pow2(log(mQ3))))*pow3(Xb))/pow3(mD32 - mQ32) - (4*dilog(1 - m32/
        mD32)*pow4(m3))/pow2(m32 - mD32) - (4*dilog(1 - m32/mQ32)*pow4(m3))/
        pow2(m32 - mQ32) + 2*pow2(Xb)*((-2*m32)/(mD32*mQ32) - (2*lMR*m32)/(
        mD32*mQ32) + (8*lMR*log(mD3))/(mD32 - mQ32) + (4*(3*m32 - 2*mD32)*log(
        mD3))/((m32 - mD32)*(mD32 - mQ32)) - (8*lMR*log(mQ3))/(mD32 - mQ32) + (
        4*(3*m32 - 2*mQ32)*log(mQ3))/((m32 - mQ32)*(-mD32 + mQ32)) + (8*(mD32 +
        mQ32)*log(mD3)*log(mQ3))/pow2(mD32 - mQ32) - (4*(3*mD32 - mQ32)*pow2(
        log(mD3)))/pow2(mD32 - mQ32) + (4*(mD32 - 3*mQ32)*pow2(log(mQ3)))/pow2(
        mD32 - mQ32) + (4*(m32 - mD32 - mQ32)*log(m3)*pow4(m3))/((m32 - mD32)*
        mD32*(m32 - mQ32)*mQ32)) + (6*log(mD3)*(-2*m32*mD32 + 2*pow4(m3) +
        pow4(mD3)))/pow2(m32 - mD32) - (4*pow2(log(mD3))*(-2*m32*mD32 + 3*pow4(
        m3) + pow4(mD3)))/pow2(m32 - mD32) + (2*(3 + 2*lMR)*log(mQ3)*(-2*m32*
        mQ32 + 2*pow4(m3) + pow4(mQ3)))/pow2(m32 - mQ32) - (4*pow2(log(mQ3))*(-
        2*m32*mQ32 + 3*pow4(m3) + pow4(mQ3)))/pow2(m32 - mQ32) + (2*(-(mD32*(
        mD32 - mQ32)*mQ32*(-2*m32 + mD32 + mQ32)*dilog(1 - m32/mD32)) + mD32*(
        mD32 - mQ32)*mQ32*(-2*m32 + mD32 + mQ32)*dilog(1 - m32/mQ32) - 2*mD32*(
        mD32 - mQ32)*mQ32*(2*m32 + 7*mD32 + 3*mQ32)*log(mD3) + lMR*(mD32 -
        mQ32)*((mD32 - mQ32)*(4*mD32*mQ32 + m32*(mD32 + mQ32)) - 4*mD32*mQ32*(
        m32 + mD32 + mQ32)*log(mD3)) + (6*mD32*mQ32 + m32*(mD32 + mQ32))*pow2(
        mD32 - mQ32) - 2*m32*(mD32 + mQ32)*log(m3)*pow2(mD32 - mQ32) + 2*mD32*
        mQ32*log(mQ3)*((mD32 - mQ32)*(2*(1 + lMR)*m32 + (3 + 2*lMR)*mD32 + (7 +
        2*lMR)*mQ32) - 4*log(mD3)*pow2(mD32 + mQ32)) + 4*mD32*mQ32*(m32*(mD32 -
        mQ32) + 2*mD32*(mD32 + mQ32))*pow2(log(mD3)) + 4*mD32*mQ32*(m32*(-mD32
        + mQ32) + 2*mQ32*(mD32 + mQ32))*pow2(log(mQ3)))*pow4(Xb))/(mD32*mQ32*
        pow4(mD32 - mQ32)) + (-6*m32*mD32*mQ32*(mD32 + mQ32) + 3*pow4(mD3)*
        pow4(mQ3) + pow4(m3)*(9*mD32*mQ32 + 2*pow4(mD3) + 2*pow4(mQ3)) - 2*(
        mD32 + mQ32)*pow6(m3))/(mD32*(-m32 + mD32)*(m32 - mQ32)*mQ32) - (2*log(
        m3)*pow4(m3)*(-4*mD32*mQ32*log(mQ3)*pow2(m32 - mD32) - 4*mD32*mQ32*log(
        mD3)*pow2(m32 - mQ32) + 2*m32*(mD32 + mQ32)*pow2(mD32 - mQ32) + pow4(
        m3)*(2*mD32*mQ32 - 4*pow4(mD3) - 4*pow4(mQ3)) + mD32*mQ32*(pow4(mD3) +
        pow4(mQ3)) + 2*lMR*mD32*mQ32*(-2*m32*(mD32 + mQ32) + 2*pow4(m3) + pow4(
        mD3) + pow4(mQ3)) + 2*(mD32 + mQ32)*pow6(m3)))/(mD32*mQ32*pow2(m32 -
        mD32)*pow2(m32 - mQ32)) + (2*lMR*(2*log(mD3)*(-2*m32*mD32 + 2*pow4(m3)
        + pow4(mD3)) + ((m32 - mD32)*(3*m32*mD32*mQ32*(mD32 + mQ32) - 3*pow4(
        mD3)*pow4(mQ3) - pow4(m3)*(3*mD32*mQ32 + pow4(mD3) + pow4(mQ3)) + (mD32
        + mQ32)*pow6(m3)))/(mD32*(m32 - mQ32)*mQ32)))/pow2(m32 - mD32));
}

double ThresholdCalculator::getDeltaLambdaYb6(int omitLogs) const
{
   using std::log;
   using himalaya::dilog;


   const double MR2 = pow2(p.scale);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mD3 = std::sqrt(mD32);
   const double mA = p.MA;
   const double beta = std::atan(p.vu/p.vd);
   const double sbeta = std::sin(beta);
   const double cbeta = std::cos(beta);
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double Yb = p.Ad(2,2) + p.mu*p.vd/p.vu;
   const double Mu = p.mu;
   const double Xb2 = pow2(Xb);
   const double mA2 = pow2(mA);
   const double Mu2 = pow2(Mu);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);
   const double lmUMR = omitLogs*log(Mu2 / MR2);
   const double eps = mQ3*0.01;

   if(std::abs(mD3 - mQ3) < eps && std::abs(mA - mQ3) < eps && std::abs(Mu - mQ3) < eps){
     const double M = std::sqrt(mQ3*mD3);
     return (30.*pow4(cbeta)*(Xb2*pow2(M)*(Xb2*(0.2 + 0.7*pow2(cbeta) - 0.1*pow2(
        sbeta)) + 0.17495808000000004*Xb*Yb*pow2(sbeta) - 0.004229546666666551*
        pow2(sbeta)*pow2(Yb)) + (Xb2*(-0.4 - 2.8*pow2(cbeta) + 
        0.5375628799999999*pow2(sbeta)) + 0.10016767999999983*Xb*Yb*pow2(sbeta)
        + 0.025041919999999957*pow2(sbeta)*pow2(Yb))*pow4(M) - 
        0.0041457066666666705*pow2(sbeta)*pow2(Yb)*pow4(Xb) + (0.8 + pow2(
        cbeta) - 0.6699863197821283*pow2(sbeta))*pow6(M) - 0.06666666666666667*
        pow2(cbeta)*pow6(Xb)))/pow6(M);
   }

   return 3*(10 - 5/pow2(cbeta) - (2*mD32)/(mQ32*pow2(cbeta)) - (4*mQ32)/(mD32*
        pow2(cbeta)) + (2*mD32)/((mD32 - Mu2)*pow2(cbeta)) + (8*Mu2)/(mD32*
        pow2(cbeta)) + (4*Mu2)/(mQ32*pow2(cbeta)) + (4*Mu2)/((mQ32 - Mu2)*pow2(
        cbeta)) + (3*pow2(sbeta))/pow2(cbeta) - (4*mA2*pow2(sbeta))/(mD32*pow2(
        cbeta)) - (2*mA2*pow2(sbeta))/(mQ32*pow2(cbeta)) + (2*pow2(Pi)*pow2(
        sbeta))/pow2(cbeta) - (4*pow2(sbeta)*pow2(Yb))/(mD32*pow2(cbeta)) - (2*
        pow2(sbeta)*pow2(Yb))/(mQ32*pow2(cbeta)) + (6*pow2(sbeta)*pow2(log(mA2/
        MR2)))/pow2(cbeta) - (96*Yb*pow2(sbeta)*pow3(Xb))/(pow2(cbeta)*pow2(
        mD32 - mQ32)) + (4*dilog(1 - Mu2/mD32)*(-2*mD32*Mu2 + pow4(mD3) - pow4(
        Mu)))/(pow2(cbeta)*pow2(mD32 - Mu2)) + (4*dilog(1 - Mu2/mQ32)*(1 - (4*
        pow4(Mu))/pow2(mQ32 - Mu2)))/pow2(cbeta) + 2*pow2(lmQ3MR)*(5 + (3 - (2*
        mD32*Mu2)/pow2(mD32 - Mu2) + 5*pow2(sbeta) + pow4(mD3)/pow2(mD32 - Mu2)
        + (-(1/pow2(mD32 - Mu2)) - 4/pow2(mQ32 - Mu2))*pow4(Mu))/pow2(cbeta)) +
        (2*lmUMR*Mu2*(-2/mQ32 + 4*((1/mQ32 + 1/(-mD32 + mQ32))/mD32 + 1/((mD32
        - Mu2)*(-mQ32 + Mu2)))*Xb2 - (2*mD32)/pow2(mD32 - Mu2) + Mu2*(3/pow2(
        mD32 - Mu2) - 6/pow2(mQ32 - Mu2) - (4*Mu2)/pow2(-(mD3*Mu2) + pow3(mD3))
        ) - (2*(-2*mD32*mQ32 + pow4(mD3) - 2*pow4(mQ3))*pow4(Xb))/(mD32*mQ32*
        pow3(mD32 - mQ32))))/pow2(cbeta) - (2*Xb2*(-12*mD32*mQ32*Mu2 + 4*mD32*
        mQ32*(mQ32 - Mu2)*dilog(1 - Mu2/mD32) + 4*mD32*mQ32*Mu2*dilog(1 - Mu2/
        mQ32) + 6*mD32*(mD32 - mQ32)*mQ32*dilog(1 - mD32/mQ32)*pow2(cbeta) + 6*
        mA2*mD32*mQ32*pow2(sbeta) + 6*mD32*mQ32*pow2(sbeta)*pow2(Yb) + 8*mQ32*
        pow4(mD3) + 4*Mu2*pow4(mD3) - 4*mQ32*pow2(cbeta)*pow4(mD3) - 2*mA2*
        pow2(sbeta)*pow4(mD3) - 4*mQ32*pow2(sbeta)*pow4(mD3) - 2*pow2(sbeta)*
        pow2(Yb)*pow4(mD3) - 2*mD32*pow4(mQ3) + 8*Mu2*pow4(mQ3) - 4*mD32*dilog(
        1 - Mu2/mQ32)*pow4(mQ3) + mD32*pow2(cbeta)*pow4(mQ3) - 4*mA2*pow2(
        sbeta)*pow4(mQ3) + 4*mD32*pow2(sbeta)*pow4(mQ3) - 4*pow2(sbeta)*pow2(
        Yb)*pow4(mQ3) - 2*pow6(mD3) + pow2(cbeta)*pow6(mD3) - 4*pow6(mQ3) + 2*
        pow2(cbeta)*pow6(mQ3)))/(mD32*mQ32*pow2(cbeta)*pow2(mD32 - mQ32)) - (2*
        (6*mD32*(mD32 - mQ32)*mQ32*dilog(1 - mD32/mQ32) - 6*mQ32*pow4(mD3) +
        15*mD32*pow4(mQ3) + pow6(mD3) + 2*pow6(mQ3))*pow6(Xb))/(mD32*mQ32*pow4(
        mD32 - mQ32)) + 2*pow2(log(mD32/mQ32))*((4*Xb*Yb*pow2(sbeta))/((mD32 -
        mQ32)*pow2(cbeta)) - (4*(-mA2 + mD32 + mQ32)*Yb*pow2(sbeta)*pow3(Xb))/(
        pow2(cbeta)*pow3(mD32 - mQ32)) + (2*Xb2*(mQ32*(-Mu2 + mQ32*(1 + 3*pow2(
        cbeta) + pow2(sbeta)) - pow2(sbeta)*pow2(Yb)) + mD32*(Mu2 - mQ32*(9 +
        3*pow2(cbeta) + 2*pow2(sbeta)) + pow2(sbeta)*pow2(Yb)) + (8 + 3*pow2(
        cbeta) + pow2(sbeta))*pow4(mD3)))/(pow2(cbeta)*pow3(mD32 - mQ32)) + (-
        2*mD32*Mu2*(1 + 2*pow2(cbeta) + 2*pow2(sbeta)) + (1 + 2*pow2(cbeta) +
        2*pow2(sbeta))*pow4(mD3) + (-1 + 2*pow2(cbeta) + 2*pow2(sbeta))*pow4(
        Mu))/(pow2(cbeta)*pow2(mD32 - Mu2)) - ((-4*mA2*pow2(sbeta)*pow2(Yb) +
        2*mQ32*(3*Mu2 + pow2(sbeta)*pow2(Yb)) + mD32*(mQ32*(14 - 8*pow2(cbeta))
        + 6*pow2(sbeta)*pow2(Yb)) + (17 + 2*pow2(cbeta))*pow4(mD3) - 2*(1 + 3*
        pow2(cbeta))*pow4(mQ3) - 3*pow4(Mu))*pow4(Xb))/(pow2(cbeta)*pow4(mD32 -
        mQ32)) - (2*(-2*mQ32*pow4(mD3) + 4*mD32*pow4(mQ3) + 3*pow6(mD3) + pow6(
        mQ3))*pow6(Xb))/pow6(mD32 - mQ32)) + (2*phixyz(mA2,mD32,mD32)*pow2(
        sbeta)*((mA2*(mA2 - 6*mD32)*Xb2*(mD32 - mQ32 + pow2(Yb)))/(pow2(mD32 -
        mQ32)*pow4(mD3)) + (2*mA2*(mA2 - 6*mD32)*Yb*pow3(Xb))/(pow2(mD32 -
        mQ32)*pow4(mD3)) + (-6*mA2*mD32 + pow4(mA))/pow4(mD3) + (mA2*(mA2 - 6*
        mD32)*pow2(Yb)*pow4(Xb))/(pow3(mD32 - mQ32)*pow4(mD3)) + (2*Xb*Yb*(-6*
        mA2*mD32 + pow4(mA)))/(-(mQ32*pow4(mD3)) + pow6(mD3)) - (deltaxyz(mA2,
        mD32,mD32)*((mQ32 - Xb2)*pow2(-(mQ3*Xb*Yb) + pow3(mQ3)) + pow4(mD3)*(
        Xb2*Yb*(6*Xb + Yb) - 3*mQ32*Xb*(Xb + 2*Yb) + 6*pow4(mQ3)) + (-4*mQ32 +
        Xb*(Xb + 2*Yb))*pow6(mD3) + mD32*(-2*mQ32*Xb2*Yb*(4*Xb + Yb) + 3*Xb*(Xb
        + 2*Yb)*pow4(mQ3) + 5*pow2(Yb)*pow4(Xb) - 4*pow6(mQ3)) + pow8(mD3)))/
        pow4(-(mD3*mQ32) + pow3(mD3))))/pow2(cbeta) + log(mD32/mQ32)*(-8 + (24*
        Xb*Yb*pow2(sbeta))/((-mD32 + mQ32)*pow2(cbeta)) + (24*(3*mD32 + mQ32)*
        Yb*pow2(sbeta)*pow3(Xb))/(pow2(cbeta)*pow3(mD32 - mQ32)) - (2*Xb*log(
        mA2/MR2)*pow2(sbeta)*(mA2*Xb*(-6*mQ32*Xb*(Xb + 2*Yb) - 2*mD32*(mQ32 -
        6*Xb*Yb) - 3*Xb2*pow2(Yb) + pow4(mD3) + pow4(mQ3)) - (mD32 - mQ32)*(
        mQ32*Xb*Yb*(2*Xb + Yb) - mD32*(Xb*Yb*(2*Xb + Yb) + 2*mQ32*(Xb + 2*Yb))
        + 3*pow2(Yb)*pow3(Xb) + (Xb + 2*Yb)*pow4(mD3) + (Xb + 2*Yb)*pow4(mQ3)))
        )/(pow2(cbeta)*pow4(mD32 - mQ32)) + (mD32*(2/mQ32 + (4*Mu2)/pow2(mD32 -
        Mu2)) - ((mA2 + 3*mQ32)*pow2(sbeta)*pow2(Yb))/(mD32*mQ32) + pow2(sbeta)
        *(-8 + (3*pow2(Yb))/mQ32) - (2*pow4(mD3))/pow2(mD32 - Mu2) + (4*pow4(
        Mu))/pow2(mD32 - Mu2))/pow2(cbeta) + (4*lmUMR*pow4(Mu)*(1/pow2(mD32 -
        Mu2) - (3*pow4(Xb))/pow4(mD32 - mQ32)))/pow2(cbeta) + (pow2(sbeta)*
        pow2(Yb)*((-5*mD32 + mQ32)*pow4(mA) + 7*mQ32*pow4(mD3) + 5*mA2*(-2*
        mD32*mQ32 + pow4(mD3) - pow4(mQ3)) - 9*mD32*pow4(mQ3) + pow6(mA) -
        pow6(mD3) - (pow4(Xb)*(-((3*mD32 + 7*mQ32)*pow4(mA)) - 9*mQ32*pow4(mD3)
        + 17*mD32*pow4(mQ3) + mA2*(10*mD32*mQ32 + 5*pow4(mD3) + 11*pow4(mQ3)) +
        pow6(mA) - 3*pow6(mD3) - 5*pow6(mQ3)))/pow2(mD32 - mQ32) + 3*pow6(mQ3)
        + (2*Xb2*((mD32 + 4*mQ32)*pow4(mA) + 8*mQ32*pow4(mD3) - 11*mD32*pow4(
        mQ3) - 2*mA2*(7*mD32*mQ32 + 4*pow4(mQ3)) - pow6(mD3) + 4*pow6(mQ3)))/(
        mD32 - mQ32)))/(mD32*mQ32*deltaxyz(mA2,mQ32,mD32)*pow2(cbeta)) - (2*(-(
        mQ32*pow4(mD3)) + 11*mD32*pow4(mQ3) + pow6(mD3) + 13*pow6(mQ3))*pow6(
        Xb))/(mQ32*pow5(-mD32 + mQ32)) + 2*lmQ3MR*((10*Xb*Yb*pow2(sbeta))/((
        mD32 - mQ32)*pow2(cbeta)) + (2*(6*mA2 - 5*mD32 - 7*mQ32)*Yb*pow2(sbeta)
        *pow3(Xb))/(pow2(cbeta)*pow3(mD32 - mQ32)) + (Xb2*(mQ32*(-4*Mu2 + mQ32*
        (22 + 16*pow2(cbeta) + 3*pow2(sbeta)) + pow2(sbeta)*(mA2 + pow2(Yb))) -
        mD32*(-4*Mu2 + mQ32*(38 + 22*pow2(cbeta) + 8*pow2(sbeta)) + pow2(sbeta)
        *(mA2 + pow2(Yb))) + (16 + 18*pow2(cbeta) + 5*pow2(sbeta))*pow4(mD3)))/
        (pow2(cbeta)*pow3(mD32 - mQ32)) + (4*(-2*mD32*Mu2*(1 + pow2(cbeta) +
        pow2(sbeta)) + (1 + pow2(cbeta) + pow2(sbeta))*pow4(mD3) + (pow2(cbeta)
        + pow2(sbeta))*pow4(Mu)))/(pow2(cbeta)*pow2(mD32 - Mu2)) - ((3*mA2*
        pow2(sbeta)*pow2(Yb) + 3*mQ32*(4*Mu2 - 3*pow2(sbeta)*pow2(Yb)) + mD32*(
        -2*mQ32*(3 + 2*pow2(cbeta)) + 3*pow2(sbeta)*pow2(Yb)) + 2*(9 + 7*pow2(
        cbeta))*pow4(mD3) - 2*(12 + 5*pow2(cbeta))*pow4(mQ3) - 6*pow4(Mu))*
        pow4(Xb))/(pow2(cbeta)*pow4(mD32 - mQ32)) + (6*mQ32*pow6(Xb))/pow4(mD32
        - mQ32)) + (2*Xb2*(-4*Mu2*pow2(sbeta)*pow2(Yb)*pow4(mQ3) + pow4(mD3)*(
        mQ32*(4*Mu2*(1 + 6*pow2(cbeta) + 2*pow2(sbeta)) + pow2(sbeta)*(2*mA2 -
        pow2(Yb))) + Mu2*pow2(sbeta)*pow2(Yb) + (22 + 13*pow2(cbeta) + 4*pow2(
        sbeta))*pow4(mQ3)) + mD32*mQ32*(Mu2*pow2(sbeta)*(-2*mA2 + pow2(Yb)) -
        mQ32*(Mu2*(18 + 13*pow2(cbeta) + 4*pow2(sbeta)) - 4*pow2(sbeta)*pow2(
        Yb)) + 4*pow4(Mu)) - (Mu2*(-2 + pow2(cbeta)) + 4*mQ32*(3 + 6*pow2(
        cbeta) + 2*pow2(sbeta)) + pow2(sbeta)*pow2(Yb))*pow6(mD3) + (-2 + pow2(
        cbeta))*pow8(mD3)))/(mD32*mQ32*(mD32 - Mu2)*pow2(cbeta)*pow2(mD32 -
        mQ32)) - (pow4(Xb)*((-mA2 + 5*mQ32)*pow2(sbeta)*pow2(Yb)*pow4(mQ3) +
        pow4(mD3)*(-(mA2*pow2(sbeta)*pow2(Yb)) - 3*mQ32*(4*Mu2 - 5*pow2(sbeta)*
        pow2(Yb)) + (90 + 76*pow2(cbeta))*pow4(mQ3)) + (-20*mQ32*(5 + 4*pow2(
        cbeta)) + pow2(sbeta)*pow2(Yb))*pow6(mD3) + mD32*(2*mA2*mQ32*pow2(
        sbeta)*pow2(Yb) + 3*(-8*Mu2 + pow2(sbeta)*(4*mA2 + pow2(Yb)))*pow4(mQ3)
        + 48*pow6(mQ3)) + (-2 + 4*pow2(cbeta))*pow8(mD3)))/(mD32*mQ32*pow2(
        cbeta)*pow4(mD32 - mQ32))) + (log(mA2/MR2)*pow2(sbeta)*(4*mD32*mQ32*
        pow2(Yb) + mA2*(mD32 + 2*mQ32)*(2*mD32 + pow2(Yb)) + (14*mQ32 - pow2(
        Yb))*pow4(mD3) - 2*pow2(Yb)*pow4(mQ3) - (2*Xb2*(2*mA2*(-2*mD32*mQ32 -
        mQ32*pow2(Yb) + pow4(mD3)) + pow2(Yb)*(-4*mD32*mQ32 + pow4(mD3) + 2*
        pow4(mQ3))))/(mD32 - mQ32) + (pow4(Xb)*(mA2*(-((10*mQ32 + pow2(Yb))*
        pow4(mD3)) + mD32*(3*mQ32*pow2(Yb) - 4*pow4(mQ3)) - 2*pow2(Yb)*pow4(
        mQ3) + 2*pow6(mD3)) + pow2(Yb)*(mQ32*pow4(mD3) - 6*mD32*pow4(mQ3) + 3*
        pow6(mD3) + 2*pow6(mQ3))))/pow3(mD32 - mQ32) + (pow2(Yb)*(-(pow6(mA)*(
        2*pow2(-(mQ3*Xb2) + pow3(mQ3)) - mD32*(-4*mQ32*Xb2 + 3*pow4(mQ3) +
        pow4(Xb)) + pow6(mD3))) + pow2(mD32 - mQ32)*(2*pow2(mQ32 - Xb2)*pow4(
        mQ3) + pow4(mD3)*(-10*mQ32*Xb2 + 11*pow4(mQ3) - 3*pow4(Xb)) - 4*mD32*
        mQ32*(-3*mQ32*Xb2 + 2*pow4(mQ3) + pow4(Xb)) + (-6*mQ32 + 2*Xb2)*pow6(
        mD3) + pow8(mD3)) + pow4(mA)*(6*pow2(mQ32 - Xb2)*pow4(mQ3) - pow4(mD3)*
        (-10*mQ32*Xb2 + pow4(mQ3) + 3*pow4(Xb)) - 2*(2*mQ32 + Xb2)*pow6(mD3) +
        mD32*(4*Xb2*pow4(mQ3) + 2*mQ32*pow4(Xb) - 6*pow6(mQ3)) + 5*pow8(mD3)) -
        mA2*(mD32 - mQ32)*(-6*pow2(mQ32 - Xb2)*pow4(mQ3) + pow4(mD3)*(-4*mQ32*
        Xb2 + 5*pow4(mQ3) - 5*pow4(Xb)) - 13*mQ32*pow6(mD3) + mD32*(-8*Xb2*
        pow4(mQ3) - mQ32*pow4(Xb) + 9*pow6(mQ3)) + 5*pow8(mD3))))/(deltaxyz(
        mA2,mQ32,mD32)*pow2(mD32 - mQ32))))/(mQ32*pow2(cbeta)*pow4(mD3)) + (
        phixyz(mA2,mQ32,mD32)*pow2(sbeta)*((4*mD32*Xb*Yb*(-3*mD32*mQ32 - mA2*(
        3*mD32 + 2*mQ32) + pow4(mA) + 2*pow4(mD3) + pow4(mQ3)))/(mD32 - mQ32) +
        (4*mD32*Yb*pow3(Xb)*(-3*mD32*mQ32 - mA2*(3*mD32 + 2*mQ32) + pow4(mA) +
        2*pow4(mD3) + pow4(mQ3)))/pow2(mD32 - mQ32) - pow2(Yb)*(-12*mD32*mQ32 -
        12*mA2*(mD32 + mQ32) + 6*pow4(mA) + 5*pow4(mD3) + 6*pow4(mQ3)) - (2*
        Xb2*pow2(Yb)*(3*(3*mD32 - 2*mQ32)*pow2(mD32 - mQ32) + (9*mD32 - 6*mQ32)
        *pow4(mA) + mA2*(-6*mD32*mQ32 - 23*pow4(mD3) + 12*pow4(mQ3))))/pow2(
        mD32 - mQ32) + 2*deltaxyz(mA2,mQ32,mD32)*((-2*mD32*Xb*Yb)/(mD32 - mQ32)
        + 2*pow2(Yb) - (2*mD32*(3*mD32 - mQ32)*Yb*pow3(Xb))/pow3(mD32 - mQ32) +
        (Xb2*(7*mD32*pow2(Yb) - 4*mQ32*pow2(Yb) + pow4(mD3)))/pow2(mD32 - mQ32)
        + (pow2(Yb)*(-9*mD32*mQ32 + 18*pow4(mD3) + 2*pow4(mQ3))*pow4(Xb))/pow4(
        mD32 - mQ32)) - (pow2(Yb)*pow4(Xb)*(2*(8*mD32 - 3*mQ32)*pow4(mA) - 49*
        mQ32*pow4(mD3) - 4*mA2*(5*mD32*mQ32 + 9*pow4(mD3) - 3*pow4(mQ3)) + 28*
        mD32*pow4(mQ3) + 27*pow6(mD3) - 6*pow6(mQ3)))/pow3(mD32 - mQ32) + (
        pow2(Yb)*(-8*(mD32 + mQ32)*pow2(mD32 - mQ32 + Xb2)*pow6(mA) + 2*pow2(
        mD32 - mQ32 + Xb2)*pow8(mA) + pow2(mD32 - mQ32)*(2*pow2(mQ32 - Xb2)*
        pow4(mQ3) + pow4(mD3)*(-10*mQ32*Xb2 + 11*pow4(mQ3) - 3*pow4(Xb)) - 4*
        mD32*mQ32*(-3*mQ32*Xb2 + 2*pow4(mQ3) + pow4(Xb)) + (-6*mQ32 + 2*Xb2)*
        pow6(mD3) + pow8(mD3)) - 2*mA2*(mD32 - mQ32)*(-4*pow2(mQ32 - Xb2)*pow4(
        mQ3) - pow4(mD3)*(2*mQ32*Xb2 + pow4(mQ3) - 3*pow4(Xb)) + (-6*mQ32 + 2*
        Xb2)*pow6(mD3) + 8*mD32*(-(Xb2*pow4(mQ3)) + pow6(mQ3)) + 3*pow8(mD3)) +
        pow4(mA)*(12*pow2(mQ32 - Xb2)*pow4(mQ3) + 8*mD32*mQ32*(mQ32*Xb2 - 2*
        pow4(mQ3) + pow4(Xb)) + pow4(mD3)*(-6*mQ32*Xb2 + 3*pow4(mQ3) + 11*pow4(
        Xb)) + (-6*mQ32 + 22*Xb2)*pow6(mD3) + 7*pow8(mD3))))/(deltaxyz(mA2,
        mQ32,mD32)*pow2(mD32 - mQ32))))/(pow2(cbeta)*pow6(mD3)) + (2*pow4(Xb)*(
        6*mD32*mQ32*dilog(1 - mD32/mQ32)*pow2(cbeta)*pow2(mD32 - mQ32) - 18*
        mQ32*Mu2*pow4(mD3) + 6*mA2*mQ32*pow2(sbeta)*pow4(mD3) + 12*mQ32*pow2(
        sbeta)*pow2(Yb)*pow4(mD3) + 12*mD32*Mu2*pow4(mQ3) + 12*mD32*Mu2*dilog(1
        - Mu2/mQ32)*pow4(mQ3) - 3*mA2*mD32*pow2(sbeta)*pow4(mQ3) - 9*mD32*pow2(
        sbeta)*pow2(Yb)*pow4(mQ3) + 85*pow4(mD3)*pow4(mQ3) - 4*dilog(1 - Mu2/
        mQ32)*pow4(mD3)*pow4(mQ3) + 114*pow2(cbeta)*pow4(mD3)*pow4(mQ3) - 2*
        mD32*mQ32*dilog(1 - Mu2/mD32)*(-2*mD32*mQ32 + 6*mQ32*Mu2 + pow4(mD3) -
        2*pow4(mQ3) - 3*pow4(Mu)) - 6*mD32*mQ32*dilog(1 - Mu2/mQ32)*pow4(Mu) -
        33*mQ32*pow6(mD3) + 2*Mu2*pow6(mD3) + 2*mQ32*dilog(1 - Mu2/mQ32)*pow6(
        mD3) - 58*mQ32*pow2(cbeta)*pow6(mD3) - mA2*pow2(sbeta)*pow6(mD3) -
        pow2(sbeta)*pow2(Yb)*pow6(mD3) - 49*mD32*pow6(mQ3) + 4*Mu2*pow6(mQ3) -
        4*mD32*dilog(1 - Mu2/mQ32)*pow6(mQ3) - 62*mD32*pow2(cbeta)*pow6(mQ3) -
        2*mA2*pow2(sbeta)*pow6(mQ3) - 2*pow2(sbeta)*pow2(Yb)*pow6(mQ3) - pow8(
        mD3) + 2*pow2(cbeta)*pow8(mD3) - 2*pow8(mQ3) + 4*pow2(cbeta)*pow8(mQ3))
        )/(mD32*mQ32*pow2(cbeta)*pow4(mD32 - mQ32)) + (2*phixyz(mA2,mQ32,mQ32)*
        pow2(sbeta)*((mA2*(2*mA2 - 11*mQ32)*Xb2*(-mD32 + mQ32 + pow2(Yb)))/(
        pow2(mD32 - mQ32)*pow4(mQ3)) + (2*mA2*(2*mA2 - 11*mQ32)*Yb*pow3(Xb))/(
        pow2(mD32 - mQ32)*pow4(mQ3)) + (-11*mA2*mQ32 + 2*pow4(mA))/pow4(mQ3) +
        (mA2*(2*mA2 - 11*mQ32)*pow2(Yb)*pow4(Xb))/(pow3(-mD32 + mQ32)*pow4(mQ3)
        ) + (2*Xb*Yb*(-11*mA2*mQ32 + 2*pow4(mA)))/(-(mD32*pow4(mQ3)) + pow6(
        mQ3)) - (deltaxyz(mA2,mQ32,mQ32)*(2*Xb2*Yb*(6*Xb + Yb)*pow4(mQ3) +
        pow4(mD3)*(2*Xb2*Yb*(2*Xb + Yb) + mQ32*Xb*(7*Xb + 12*Yb) + 12*pow4(mQ3)
        ) + 9*mQ32*pow2(Yb)*pow4(Xb) - 2*(4*mQ32 + Xb*(Xb + 2*Yb))*pow6(mD3) +
        Xb*(3*Xb + 4*Yb)*pow6(mQ3) - 2*mD32*(2*mQ32*Xb2*Yb*(4*Xb + Yb) + 2*Xb*(
        2*Xb + 3*Yb)*pow4(mQ3) + pow2(Yb)*pow4(Xb) + 4*pow6(mQ3)) + 2*pow8(mD3)
        + 2*pow8(mQ3)))/pow4(-(mD32*mQ3) + pow3(mQ3))))/pow2(cbeta) + lmQ3MR*(-
        20 + (48*Yb*pow2(sbeta)*pow3(Xb))/(pow2(cbeta)*pow2(mD32 - mQ32)) + (4*
        lmUMR*(1/pow2(mD32 - Mu2) + 2/pow2(mQ32 - Mu2))*pow4(Mu))/pow2(cbeta) +
        (mD32*(2/mQ32 + (4*Mu2)/pow2(mD32 - Mu2)) + (4*mQ32*Mu2)/pow2(mQ32 -
        Mu2) - 20*pow2(sbeta) + (3*pow2(sbeta)*pow2(Yb))/mQ32 + (2*(-mA2 +
        mQ32)*pow2(sbeta)*pow2(Yb))/pow4(mD3) - (2*pow4(mD3))/pow2(mD32 - Mu2)
        + (2*pow4(mQ3))/pow2(mQ32 - Mu2) + (-(mA2*pow2(sbeta)*pow2(Yb)) + 4*
        pow4(mQ3))/(mD32*mQ32) + (4/pow2(mD32 - Mu2) + 6/pow2(mQ32 - Mu2))*
        pow4(Mu))/pow2(cbeta) + (2*(-5*mD32*mQ32 + pow4(mD3) - 2*pow4(mQ3))*
        pow6(Xb))/(mD32*mQ32*pow3(mD32 - mQ32)) - (pow4(Xb)*(2*(-mA2 + mQ32)*
        pow2(sbeta)*pow2(Yb)*pow4(mQ3) + pow4(mD3)*(-(mA2*pow2(sbeta)*pow2(Yb))
        + mQ32*(-12*Mu2 + 11*pow2(sbeta)*pow2(Yb)) + (86 + 68*pow2(cbeta))*
        pow4(mQ3)) + (-64*mQ32*(1 + pow2(cbeta)) + pow2(sbeta)*pow2(Yb))*pow6(
        mD3) + mD32*(3*mA2*mQ32*pow2(sbeta)*pow2(Yb) - 2*pow2(sbeta)*pow2(Yb)*
        pow4(mQ3) + (4 - 8*pow2(cbeta))*pow6(mQ3)) + (-2 + 4*pow2(cbeta))*pow8(
        mD3)))/(mQ32*pow2(cbeta)*pow3(mD32 - mQ32)*pow4(mD3)) + (pow2(sbeta)*
        pow2(Yb)*(pow6(mA)*(2*pow2(-(mQ3*Xb2) + pow3(mQ3)) - mD32*(-4*mQ32*Xb2
        + 3*pow4(mQ3) + pow4(Xb)) + pow6(mD3)) + pow4(mA)*(-6*pow2(mQ32 - Xb2)*
        pow4(mQ3) + pow4(mD3)*(-10*mQ32*Xb2 + pow4(mQ3) + 3*pow4(Xb)) + 2*(2*
        mQ32 + Xb2)*pow6(mD3) + mD32*(-4*Xb2*pow4(mQ3) - 2*mQ32*pow4(Xb) + 6*
        pow6(mQ3)) - 5*pow8(mD3)) - pow2(mD32 - mQ32)*(2*pow2(mQ32 - Xb2)*pow4(
        mQ3) + pow4(mD3)*(-10*mQ32*Xb2 + 11*pow4(mQ3) - 3*pow4(Xb)) - 4*mD32*
        mQ32*(-3*mQ32*Xb2 + 2*pow4(mQ3) + pow4(Xb)) + (-6*mQ32 + 2*Xb2)*pow6(
        mD3) + pow8(mD3)) + mA2*(mD32 - mQ32)*(-6*pow2(mQ32 - Xb2)*pow4(mQ3) +
        pow4(mD3)*(-4*mQ32*Xb2 + 5*pow4(mQ3) - 5*pow4(Xb)) - 13*mQ32*pow6(mD3)
        + mD32*(-8*Xb2*pow4(mQ3) - mQ32*pow4(Xb) + 9*pow6(mQ3)) + 5*pow8(mD3)))
        )/(mQ32*deltaxyz(mA2,mQ32,mD32)*pow2(cbeta)*pow2(mD32 - mQ32)*pow4(mD3)
        ) + (2*Xb2*(2*(-mA2 + mQ32)*(mQ32 - Mu2)*Mu2*pow2(sbeta)*pow2(Yb)*pow4(
        mQ3) - 2*mD32*mQ32*(mQ32 - Mu2)*(-(mA2*Mu2*pow2(sbeta)*pow2(Yb)) +
        mQ32*(-mA2 + Mu2)*pow2(sbeta)*pow2(Yb) + (Mu2*(-2 + pow2(cbeta)) +
        pow2(sbeta)*pow2(Yb))*pow4(mQ3)) + pow6(mD3)*((-14*Mu2 + pow2(sbeta)*
        pow2(Yb))*pow4(mQ3) - 2*mQ32*(-3 + pow2(sbeta))*pow4(Mu) - pow2(sbeta)*
        pow2(Yb)*pow4(Mu) + pow2(cbeta)*(15*Mu2*pow4(mQ3) - 8*mQ32*pow4(Mu) -
        7*pow6(mQ3)) + 2*pow2(sbeta)*pow6(mQ3)) + (mQ32*(7*Mu2*pow2(cbeta) +
        pow2(sbeta)*(2*Mu2 - pow2(Yb))) + Mu2*pow2(sbeta)*pow2(Yb) - 2*(-3 + 4*
        pow2(cbeta) + pow2(sbeta))*pow4(mQ3) + (-2 + pow2(cbeta))*pow4(Mu))*
        pow8(mD3) + pow4(mD3)*(mQ32*Mu2*(2*mA2 + Mu2)*pow2(sbeta)*pow2(Yb) -
        pow4(mQ3)*(pow2(sbeta)*(2*mA2*pow2(Yb) + 3*Mu2*pow2(Yb) - 2*pow4(Mu)) +
        7*pow2(cbeta)*pow4(Mu)) + (Mu2*(8 + 5*pow2(cbeta) - 2*pow2(sbeta)) + 2*
        pow2(sbeta)*pow2(Yb))*pow6(mQ3) + 2*(-2 + pow2(cbeta))*pow8(mQ3)) + (
        mQ32 - Mu2)*(-2 + pow2(cbeta))*power10(mD3)))/(mQ32*(mD32 - Mu2)*(mQ32
        - Mu2)*pow2(cbeta)*pow2(mD32 - mQ32)*pow4(mD3))));
}

/**
 * 2-loop threshold correction to lambda O(at^2) from SUSYHD
 * [1504.05200] Eq.(21).
 *
 * @note Only valid for mQ3 = mU3!
 *
 * @note Multiplied by a factor 2 in order to match the convention of
 * LTM.
 *
 * @param omitLogs factor which multiplies the log(mst^2/Q^2) terms
 *
 * @return 2-loop threshold correction O(at^2)
 */
double ThresholdCalculator::getDeltaLambdaYt6_SUSYHD(int omitLogs) const
{
   using std::log;

   const double mQ32 = p.mq2(2,2);
   const double mU32 = p.mu2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mU3 = std::sqrt(mU32);
   const double Xt = p.Au(2,2) - p.mu*p.vd/p.vu;
   const double mst2 = mQ3 * mU3;
   const double mst = std::sqrt(mst2);
   const double xt = Xt/mst;
   const double muhat = p.mu/mst;
   const double tb = p.vu/p.vd;
   const double beta = std::atan(tb);
   const double cb2 = pow2(std::cos(beta));
   const double sb2 = pow2(std::sin(beta));
   const double s2b = std::sin(2*beta);
   const double Q2 = pow2(p.scale);
   const double lmst2Q2 = omitLogs*std::log(mst2/Q2);

   const double dlam =
      (3*(0.5 - 8.34993159891064*cb2 - 4*lmst2Q2 + 13*cb2*lmst2Q2 +
       12.*cb2*(0.04173653333333327 + lmst2Q2)*xt*((2*muhat)/s2b + xt) -
       8*f1HD(muhat) + 4*f3HD(muhat) + 3*pow2(lmst2Q2) -
       3*cb2*pow2(lmst2Q2) + 6*pow2(muhat) - 6*lmst2Q2*pow2(muhat) -
       2*f1HD(muhat)*pow2(muhat) + 3*f2HD(muhat)*pow2(muhat) +
       3.*cb2*(0.04173653333333327 + lmst2Q2)*pow2((2*muhat)/s2b + xt) +
       pow2(xt)*(-7 + 19.6878144*cb2 + 27*lmst2Q2 - 24*cb2*lmst2Q2 +
          4*f1HD(muhat) - 4*f2HD(muhat) - 6*pow2(muhat) +
          6*lmst2Q2*pow2(muhat) - 6*f1HD(muhat)*pow2(muhat) -
          6*f2HD(muhat)*pow2(muhat) -
          3.*cb2*(0.007049244444444251 + lmst2Q2)*pow2((2*muhat)/s2b + xt))
        - 2.*cb2*(-0.4373952000000001 + lmst2Q2)*((2*muhat)/s2b + xt)*
        pow3(xt) + ((22 - 25*cb2 - 26*lmst2Q2 + 24*cb2*lmst2Q2 -
            2*f1HD(muhat) + 2*f2HD(muhat) + 4*pow2(muhat) -
            4*lmst2Q2*pow2(muhat) + 4*f1HD(muhat)*pow2(muhat) +
            2*f2HD(muhat)*pow2(muhat) +
            (2.*cb2*(-0.04145706666666671 + lmst2Q2)*
               pow2(2.*muhat + s2b*xt))/pow2(s2b))*pow4(xt))/4. -
       ((-1 + cb2)*(-1 + lmst2Q2)*pow6(xt))/2.))/sb2;

   return 2*dlam;
}

double ThresholdCalculator::getDeltaLambdaYt6(int omitLogs) const
{
   using std::log;
   using himalaya::dilog;

   const double MR2 = pow2(p.scale);
   double mQ32 = p.mq2(2,2);
   double mU32 = p.mu2(2,2);
   const double mA = p.MA;
   const double beta = std::atan(p.vu/p.vd);
   const double sbeta = std::sin(beta);
   const double cbeta = std::cos(beta);
   const double Xt = p.Au(2,2) - p.mu*p.vd/p.vu;
   const double Yt = p.Au(2,2) + p.mu*p.vu/p.vd;
   const double Xt2 = pow2(Xt);
   const double Xt4 = pow2(Xt2);
   const double Xt6 = pow3(Xt2);
   const double mA2 = pow2(mA);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);

   if (is_equal_rel(mQ32, mU32, 0.1) &&
       (is_equal_rel(mQ32, mA2, 0.1) || is_equal_rel(mU32, mA2, 0.1))) {
      return getDeltaLambdaYt6_SUSYHD(omitLogs);
   }

   mQ32 = p.mq2(2,2) * (1.02);
   mU32 = p.mu2(2,2) * (0.98);
   const double mQ3 = std::sqrt(mQ32);
   const double mU3 = std::sqrt(mU32);
   const double Mu = p.mu * std::sqrt(1.03);
   const double Mu2 = pow2(Mu);
   const double lmUMR = omitLogs*log(Mu2 / MR2);

   return (18*pow2(cbeta)*pow2(log(mA2/MR2)))/pow2(sbeta) - (288*Yt*pow2(cbeta)*
        pow3(Xt))/(pow2(mQ32 - mU32)*pow2(sbeta)) - (6*Xt6*(11*mQ32*mU32 + 2*
        mQ32*mU32*dilog(1 - mQ32/mU32) - 4*mQ32*mU32*dilog(1 - mU32/mQ32) + 2*
        pow4(mQ3) - pow4(mU3)))/(mQ32*mU32*pow3(mQ32 - mU32)) + (3*phixyz(mA2,
        mU32,mQ32)*pow2(cbeta)*((4*(mA2 + mQ32 - mU32)*Xt*Yt)/(mQ32 - mU32) +
        pow2(Yt) + (2*mA2*Xt2*pow2(Yt))/pow2(mQ32 - mU32) - (pow2(mA2 + mQ32 -
        mU32)*pow2(-mQ32 + mU32 + Xt2)*pow2(Yt))/(deltaxyz(mA2,mU32,mQ32)*pow2(
        mQ32 - mU32)) - ((4*mA2 + 3*mQ32 - 3*mU32)*Xt4*pow2(Yt))/pow3(mQ32 -
        mU32) - (4*(mA2 + mQ32 - mU32)*Yt*pow3(Xt))/pow2(mQ32 - mU32) + (2*Xt2*
        deltaxyz(mA2,mU32,mQ32)*(-4*mU32*Xt*Yt - 2*mQ32*(mU32 - 2*Xt*Yt) + 3*
        Xt2*pow2(Yt) + pow4(mQ3) + pow4(mU3)))/pow4(mQ32 - mU32)))/(mQ32*pow2(
        sbeta)) + (6*Xt2*(12*mQ32*Mu2*mU32 + 4*mQ32*(-mQ32 + Mu2)*mU32*dilog(1
        - mQ32/Mu2) - 4*mQ32*Mu2*mU32*dilog(1 - mU32/Mu2) - 6*mA2*mQ32*mU32*
        pow2(cbeta) + 6*mQ32*mU32*(-mQ32 + mU32)*dilog(1 - mQ32/mU32)*pow2(
        sbeta) - 6*mQ32*mU32*pow2(cbeta)*pow2(Yt) - 8*Mu2*pow4(mQ3) + 2*mU32*
        pow4(mQ3) + 4*mU32*dilog(1 - mU32/Mu2)*pow4(mQ3) + 4*mA2*pow2(cbeta)*
        pow4(mQ3) - 4*mU32*pow2(cbeta)*pow4(mQ3) - mU32*pow2(sbeta)*pow4(mQ3) +
        4*pow2(cbeta)*pow2(Yt)*pow4(mQ3) - 8*mQ32*pow4(mU3) - 4*Mu2*pow4(mU3) +
        2*mA2*pow2(cbeta)*pow4(mU3) + 4*mQ32*pow2(cbeta)*pow4(mU3) + 4*mQ32*
        pow2(sbeta)*pow4(mU3) + 2*pow2(cbeta)*pow2(Yt)*pow4(mU3) + 4*pow6(mQ3)
        - 2*pow2(sbeta)*pow6(mQ3) + 2*pow6(mU3) - pow2(sbeta)*pow6(mU3)))/(
        mQ32*mU32*pow2(mQ32 - mU32)*pow2(sbeta)) + (3*log(mA2/MR2)*pow2(cbeta)*
        (-12*lmQ3MR + (pow2(-mQ32 + mU32 + Xt2)*pow2(Yt)*(-pow2(mQ32 - mU32) +
        pow4(mA)))/(deltaxyz(mA2,mU32,mQ32)*pow2(-(mQ32*mU3) + pow3(mU3))) + (
        2*Xt2*(-(mA2*(mQ32*(4*mU32 + pow2(Yt)) - 2*pow4(mU3))) + pow2(Yt)*(-3*
        mQ32*mU32 + pow4(mQ3) + pow4(mU3))))/(mQ32*(mQ32 - mU32)*pow4(mU3)) + (
        mA2*(mU32*(2*mU32 + pow2(Yt)) + mQ32*(4*mU32 + pow2(Yt))) - pow2(Yt)*
        pow4(mQ3) - pow2(Yt)*pow4(mU3) + mQ32*(3*mU32*pow2(Yt) + 14*pow4(mU3)))
        /(mQ32*pow4(mU3)) + (pow2(Yt)*(-(mA2*(3*mQ32 + 5*mU32)*pow2(mQ32 -
        mU32)) + pow4(mA)*(4*mQ32*mU32 + 3*pow4(mQ3) + 5*pow4(mU3)) + pow4(mQ32
        - mU32) - (mQ32 + mU32)*pow6(mA) - (2*Xt2*(pow4(mA)*(2*mQ32*mU32 + 3*
        pow4(mQ3) - pow4(mU3)) + pow4(mQ32 - mU32) - mQ32*pow6(mA) + mA2*(2*
        mU32*pow4(mQ3) + mQ32*pow4(mU3) - 3*pow6(mQ3))))/(mQ32 - mU32) - (Xt4*(
        -3*(mQ32 + mU32)*pow4(mA) + 3*mU32*pow4(mQ3) + mQ32*pow4(mU3) + mA2*(3*
        pow4(mQ3) + 5*pow4(mU3)) + pow6(mA) - pow6(mQ3) - 3*pow6(mU3)))/(mQ32 -
        mU32)))/(mQ32*deltaxyz(mA2,mQ32,mU32)*pow4(mU3)) + (Xt4*(4*mU32*(-3*
        mU32 + pow2(Yt))*pow4(mQ3) - pow2(Yt)*pow6(mQ3) + mA2*((4*mU32 + pow2(
        Yt))*pow4(mQ3) + pow2(Yt)*pow4(mU3) + 2*mQ32*(-(mU32*pow2(Yt)) + 5*
        pow4(mU3)) - 2*pow6(mU3)) + 12*mQ32*pow6(mU3) - 3*pow2(Yt)*pow6(mU3)))/
        (mQ32*pow3(mQ32 - mU32)*pow4(mU3))))/pow2(sbeta) + (3*(-4*dilog(1 -
        mQ32/Mu2)*(-2*mQ32*Mu2 + pow4(mQ3) - 3*pow4(Mu)) - ((mQ32 - Mu2)*(-4*
        mQ32*(mQ32 - Mu2)*mU32*dilog(1 - mU32/Mu2)*(2*Mu2*mU32 + pow4(Mu) -
        pow4(mU3)) + (Mu2 - mU32)*((mQ32 - Mu2)*(Mu2 - mU32)*pow2(cbeta)*(2*
        mA2*(2*mQ32 + mU32) + 2*mU32*pow2(Yt) + mQ32*(-(mU32*(3 + 2*pow2(Pi)))
        + 4*pow2(Yt))) + 17*Mu2*mU32*pow4(mQ3) - 10*Mu2*mU32*pow2(sbeta)*pow4(
        mQ3) - 21*mQ32*mU32*pow4(Mu) + 10*mQ32*mU32*pow2(sbeta)*pow4(Mu) - 12*
        pow4(mQ3)*pow4(Mu) + 13*mQ32*Mu2*pow4(mU3) - 10*mQ32*Mu2*pow2(sbeta)*
        pow4(mU3) - 3*pow4(mQ3)*pow4(mU3) + 10*pow2(sbeta)*pow4(mQ3)*pow4(mU3)
        - 6*pow4(Mu)*pow4(mU3) + 4*Mu2*pow6(mQ3) - 4*mU32*pow6(mQ3) + 8*mQ32*
        pow6(Mu) + 4*mU32*pow6(Mu) - 2*mQ32*pow6(mU3) + 2*Mu2*pow6(mU3))))/(
        mQ32*mU32*pow2(Mu2 - mU32))))/(pow2(mQ32 - Mu2)*pow2(sbeta)) + (6*
        phixyz(mA2,mQ32,mQ32)*pow2(cbeta)*((mA2*(mA2 - 7*mQ32)*Xt2*(mQ32 - mU32
        + pow2(Yt)))/(pow2(mQ32 - mU32)*pow4(mQ3)) + (mA2*(mA2 - 7*mQ32)*Xt4*
        pow2(Yt))/(pow3(mQ32 - mU32)*pow4(mQ3)) + (2*mA2*(mA2 - 7*mQ32)*Yt*
        pow3(Xt))/(pow2(mQ32 - mU32)*pow4(mQ3)) + (-7*mA2*mQ32 + pow4(mA))/
        pow4(mQ3) + (2*Xt*Yt*(-7*mA2*mQ32 + pow4(mA)))/(-(mU32*pow4(mQ3)) +
        pow6(mQ3)) - (deltaxyz(mA2,mQ32,mQ32)*((mU32 - Xt2)*pow2(-(mU3*Xt*Yt) +
        pow3(mU3)) + pow4(mQ3)*(Xt2*Yt*(10*Xt + Yt) - mU32*Xt*(5*Xt + 6*Yt) +
        6*pow4(mU3)) + 2*(-2*mU32 + Xt*(Xt + Yt))*pow6(mQ3) - 2*mQ32*(mU32*Xt2*
        Yt*(6*Xt + Yt) - 4*Xt4*pow2(Yt) - Xt*(2*Xt + 3*Yt)*pow4(mU3) + 2*pow6(
        mU3)) + pow8(mQ3)))/pow4(-(mQ3*mU32) + pow3(mQ3))))/pow2(sbeta) + (12*
        pow2(lmUMR)*(2*mQ32*Mu2*mU32*(mQ32 + mU32) - pow4(mQ3)*pow4(mU3) +
        pow4(Mu)*(-4*mQ32*mU32 + pow4(mU3)) - 2*mU32*pow6(Mu) + 2*pow8(Mu)))/(
        pow2(mQ32 - Mu2)*pow2(Mu2 - mU32)*pow2(sbeta)) + (6*pow2(lmQ3MR)*(5*(
        pow2(cbeta) + pow2(sbeta))*pow4(Mu)*pow4(mU3) + pow4(mQ3)*(-2*Mu2*mU32*
        (-4 + 5*pow2(cbeta) + 5*pow2(sbeta)) + (-2 + 5*pow2(cbeta) + 5*pow2(
        sbeta))*pow4(Mu) + (-4 + 5*pow2(cbeta) + 5*pow2(sbeta))*pow4(mU3)) -
        10*mU32*(pow2(cbeta) + pow2(sbeta))*pow6(Mu) - 2*mQ32*(-2*mU32*(-4 + 5*
        pow2(cbeta) + 5*pow2(sbeta))*pow4(Mu) + Mu2*(-4 + 5*pow2(cbeta) + 5*
        pow2(sbeta))*pow4(mU3) + (-2 + 5*pow2(cbeta) + 5*pow2(sbeta))*pow6(Mu))
        + (2 + 5*pow2(cbeta) + 5*pow2(sbeta))*pow8(Mu)))/(pow2(mQ32 - Mu2)*
        pow2(Mu2 - mU32)*pow2(sbeta)) + (6*lmUMR*((4*(mQ32 - Mu2)*(Mu2 - mU32)*
        Xt2*((Mu2 - mU32)*mU32*pow4(Mu) + pow4(mQ3)*(-(Mu2*mU32) + 2*pow4(Mu))
        + mQ32*(mU32*pow4(Mu) - 2*pow6(Mu))))/(mQ32 - mU32) - Mu2*(2*mU32*pow2(
        Mu2 - mU32)*pow4(Mu) + (-3*Mu2*mU32 + 4*pow4(Mu) + 2*pow4(mU3))*pow6(
        mQ3) + pow4(mQ3)*(8*mU32*pow4(Mu) - 8*Mu2*pow4(mU3) - 8*pow6(Mu) + 2*
        pow6(mU3)) + mQ32*(-2*pow4(Mu)*pow4(mU3) - mU32*pow6(Mu) + 2*Mu2*pow6(
        mU3) + 4*pow8(Mu))) - (2*Mu2*Xt4*(-(pow2(Mu2 - mU32)*pow4(Mu)*pow4(mU3)
        ) + pow6(mQ3)*(8*mU32*pow4(Mu) - 9*Mu2*pow4(mU3) - 4*pow6(Mu) + 2*pow6(
        mU3)) + (-3*Mu2*mU32 + 2*pow4(Mu) + 2*pow4(mU3))*pow8(mQ3) + mQ32*(-5*
        pow4(mU3)*pow6(Mu) + 2*pow4(Mu)*pow6(mU3) + 2*mU32*pow8(Mu)) + pow4(
        mQ3)*(7*pow4(Mu)*pow4(mU3) - 5*mU32*pow6(Mu) + 2*pow8(Mu) - pow8(mU3)))
        )/pow3(mQ32 - mU32)))/(mQ32*mU32*pow2(mQ32 - Mu2)*pow2(Mu2 - mU32)*
        pow2(sbeta)) + (6*phixyz(mA2,mU32,mU32)*pow2(cbeta)*((mA2*(mA2 - 6*
        mU32)*Xt2*(-mQ32 + mU32 + pow2(Yt)))/(pow2(mQ32 - mU32)*pow4(mU3)) + (
        mA2*(mA2 - 6*mU32)*Xt4*pow2(Yt))/(pow3(-mQ32 + mU32)*pow4(mU3)) + (2*
        mA2*(mA2 - 6*mU32)*Yt*pow3(Xt))/(pow2(mQ32 - mU32)*pow4(mU3)) + (-6*
        mA2*mU32 + pow4(mA))/pow4(mU3) + (2*Xt*Yt*(-6*mA2*mU32 + pow4(mA)))/(-(
        mQ32*pow4(mU3)) + pow6(mU3)) - (deltaxyz(mA2,mU32,mU32)*(5*mU32*Xt4*
        pow2(Yt) + Xt2*Yt*(6*Xt + Yt)*pow4(mU3) + pow4(mQ3)*(Xt2*Yt*(2*Xt + Yt)
        + 3*mU32*Xt*(Xt + 2*Yt) + 6*pow4(mU3)) - (4*mU32 + Xt*(Xt + 2*Yt))*
        pow6(mQ3) + Xt*(Xt + 2*Yt)*pow6(mU3) - mQ32*(2*mU32*Xt2*Yt*(4*Xt + Yt)
        + Xt4*pow2(Yt) + 3*Xt*(Xt + 2*Yt)*pow4(mU3) + 4*pow6(mU3)) + pow8(mQ3)
        + pow8(mU3)))/pow4(-(mQ32*mU3) + pow3(mU3))))/pow2(sbeta) + pow2(log(
        mU32/mQ32))*((-24*Xt*Yt*pow2(cbeta))/((mQ32 - mU32)*pow2(sbeta)) + (24*
        (-mA2 + mQ32 + mU32)*Yt*pow2(cbeta)*pow3(Xt))/(pow2(sbeta)*pow3(mQ32 -
        mU32)) + 6*(2 + (2*pow2(cbeta))/pow2(sbeta) + (2*Mu2*mU32)/(pow2(Mu2 -
        mU32)*pow2(sbeta)) - pow4(mU3)/(pow2(Mu2 - mU32)*pow2(sbeta))) + (Xt2*(
        12*pow2(cbeta)*pow2(Mu2 - mU32)*(-mQ32 + mU32 + pow2(Yt)) - 6*mQ32*(2*
        Mu2*mU32*(2 - 9*pow2(sbeta)) + 9*pow2(sbeta)*pow4(Mu) + (-2 + 9*pow2(
        sbeta))*pow4(mU3)) + 6*mU32*(-2*Mu2*mU32*(14 + 9*pow2(sbeta)) + (16 +
        9*pow2(sbeta))*pow4(Mu) + (14 + 9*pow2(sbeta))*pow4(mU3))))/(pow2(mQ32
        - mU32)*pow2(Mu2 - mU32)*pow2(sbeta)) - (6*(3*mQ32 + 5*mU32)*Xt6)/pow4(
        mQ32 - mU32) + (6*Xt4*(2*mQ32*pow2(Mu2 - mU32)*(mU32*(-8 + pow2(sbeta))
        - pow2(cbeta)*pow2(Yt)) + 4*mA2*pow2(cbeta)*pow2(Yt)*pow4(mU3) + pow4(
        mQ3)*(2*Mu2*mU32*(1 - 3*pow2(sbeta)) + 3*pow2(sbeta)*pow4(Mu) + (-1 +
        3*pow2(sbeta))*pow4(mU3)) - pow4(Mu)*(-4*mA2*pow2(cbeta)*pow2(Yt) + 6*
        mU32*pow2(cbeta)*pow2(Yt) + (16 + 5*pow2(sbeta))*pow4(mU3)) - 6*pow2(
        cbeta)*pow2(Yt)*pow6(mU3) + 2*Mu2*(-4*mA2*mU32*pow2(cbeta)*pow2(Yt) +
        6*pow2(cbeta)*pow2(Yt)*pow4(mU3) + 5*(3 + pow2(sbeta))*pow6(mU3)) - 15*
        pow8(mU3) - 5*pow2(sbeta)*pow8(mU3)))/(pow2(Mu2 - mU32)*pow2(sbeta)*
        pow4(mQ32 - mU32))) + (3*phixyz(mA2,mQ32,mU32)*pow2(cbeta)*pow2(Yt)*(-
        3*(-2*mA2*(mQ32 + mU32) + pow2(mQ32 - mU32) + pow4(mA)) + (2*Xt2*((3*
        mQ32 - 5*mU32)*pow2(mQ32 - mU32) + (3*mQ32 - 5*mU32)*pow4(mA) + mA2*(4*
        mQ32*mU32 - 6*pow4(mQ3) + 14*pow4(mU3))))/pow2(mQ32 - mU32) + 2*
        deltaxyz(mA2,mQ32,mU32)*(1 - (2*(mQ32 - 2*mU32)*Xt2)/pow2(mQ32 - mU32)
        + (Xt4*(-5*mQ32*mU32 + pow4(mQ3) + 12*pow4(mU3)))/pow4(mQ32 - mU32)) +
        (Xt4*(3*(mQ32 - 3*mU32)*pow4(mA) - 15*mU32*pow4(mQ3) - 6*mA2*(-2*mQ32*
        mU32 + pow4(mQ3) - 3*pow4(mU3)) + 29*mQ32*pow4(mU3) + 3*pow6(mQ3) - 17*
        pow6(mU3)))/pow3(-mQ32 + mU32) + (-4*(mQ32 + mU32)*pow2(-mQ32 + mU32 +
        Xt2)*pow6(mA) - 4*mA2*pow2(mQ32 - mU32)*(-((mU32 + 2*Xt2)*pow4(mQ3)) +
        mQ32*(Xt4 - pow4(mU3)) + mU32*(Xt4 + pow4(mU3)) + pow6(mQ3)) + pow2(-
        mQ32 + mU32 + Xt2)*pow8(mA) + pow2(mQ32 - mU32)*(-3*Xt4*pow4(mU3) - 2*
        mQ32*mU32*(3*mU32*Xt2 + Xt4 + 2*pow4(mU3)) + pow4(mQ3)*(6*mU32*Xt2 +
        Xt4 + 6*pow4(mU3)) - 2*(2*mU32 + Xt2)*pow6(mQ3) + 2*Xt2*pow6(mU3) +
        pow8(mQ3) + pow8(mU3)) + 2*pow4(mA)*((2*mU32*Xt2 + 3*Xt4)*pow4(mQ3) +
        3*Xt4*pow4(mU3) + mQ32*(2*mU32*Xt4 - 2*Xt2*pow4(mU3)) - 2*(2*mU32 + 3*
        Xt2)*pow6(mQ3) + 6*Xt2*pow6(mU3) + 3*pow8(mQ3) + pow8(mU3)))/(deltaxyz(
        mA2,mQ32,mU32)*pow2(mQ32 - mU32))))/(pow2(sbeta)*pow6(mU3)) + log(mU32/
        mQ32)*((72*Xt*Yt*pow2(cbeta))/((mQ32 - mU32)*pow2(sbeta)) - (6*(mA2 +
        mQ32 - mU32)*pow2(cbeta)*pow2(-mQ32 + mU32 + Xt2)*pow2(Yt))/(deltaxyz(
        mA2,mU32,mQ32)*pow2(mQ32 - mU32)*pow2(sbeta)) - (72*(mQ32 + 3*mU32)*Yt*
        pow2(cbeta)*pow3(Xt))/(pow2(sbeta)*pow3(mQ32 - mU32)) + (6*log(mA2/MR2)
        *pow2(cbeta)*(-3 - (2*Xt*Yt)/(mQ32 - mU32) - (Xt2*(mA2 - 5*mQ32 + 5*
        mU32 + pow2(Yt)))/pow2(mQ32 - mU32) + (2*(6*mA2 - mQ32 + mU32)*Yt*pow3(
        Xt))/pow3(mQ32 - mU32) + (3*Xt4*(mA2*(2*mQ32 + pow2(Yt)) - (mQ32 -
        mU32)*(mQ32 + mU32 + pow2(Yt))))/pow4(mQ32 - mU32)))/pow2(sbeta) + (6*
        Xt6*(-3*mQ32*mU32 - 10*pow4(mQ3) + pow4(mU3)))/(mQ32*pow4(mQ32 - mU32))
        + (3*pow2(cbeta)*pow2(Yt)*(pow3(mQ32 - mU32) - (mQ32 + 5*mU32)*pow4(mA)
        - (2*Xt2*((2*mQ32 - mU32)*pow2(mQ32 - mU32) + (2*mQ32 + mU32)*pow4(mA)
        - 4*mA2*(2*mQ32*mU32 + pow4(mQ3))))/(mQ32 - mU32) - mA2*(4*mQ32*mU32 +
        pow4(mQ3) - 5*pow4(mU3)) + pow6(mA) - (Xt4*(-((5*mQ32 + 3*mU32)*pow4(
        mA)) + 11*mU32*pow4(mQ3) - 5*mQ32*pow4(mU3) + mA2*(4*mQ32*mU32 + 7*
        pow4(mQ3) + 5*pow4(mU3)) + pow6(mA) - 3*pow6(mQ3) - 3*pow6(mU3)))/pow2(
        mQ32 - mU32)))/(mQ32*mU32*deltaxyz(mA2,mQ32,mU32)*pow2(sbeta)) + (6*
        Xt2*((pow2(cbeta)*(2*mA2*mQ32*mU32 + (mU32 + 2*pow2(Yt))*pow4(mQ3) +
        mQ32*(mU32*pow2(Yt) - 5*pow4(mU3)) - pow2(Yt)*pow4(mU3)))/mU32 - (Mu2*(
        Mu2 - mU32)*(-2 + pow2(sbeta))*pow4(mU3) + pow4(mQ3)*(Mu2*mU32*(-16 +
        5*pow2(sbeta)) + (13 + 16*pow2(sbeta))*pow4(Mu) - 3*(3 + 7*pow2(sbeta))
        *pow4(mU3)) + (-(Mu2*(13 + 16*pow2(sbeta))) + mU32*(19 + 16*pow2(sbeta)
        ))*pow6(mQ3) + mQ32*(mU32*(5 - 21*pow2(sbeta))*pow4(Mu) + Mu2*(7 + 20*
        pow2(sbeta))*pow4(mU3) - 4*pow6(Mu) + (-2 + pow2(sbeta))*pow6(mU3)))/((
        mQ32 - Mu2)*(Mu2 - mU32))))/(mQ32*pow2(mQ32 - mU32)*pow2(sbeta)) - (3*(
        pow2(-(Mu*mU32) + pow3(Mu))*(-(mA2*pow2(cbeta)*pow2(Yt)) + 3*mU32*pow2(
        cbeta)*pow2(Yt) + 2*pow4(mU3)) + pow4(mQ3)*((mU32*(-9 + 5*pow2(cbeta) +
        8*pow2(sbeta)) + pow2(cbeta)*pow2(Yt))*pow4(Mu) + pow2(cbeta)*pow2(Yt)*
        pow4(mU3) - 2*Mu2*(mU32*pow2(cbeta)*pow2(Yt) + (-2 + 5*pow2(cbeta) + 8*
        pow2(sbeta))*pow4(mU3)) + (-1 + 5*pow2(cbeta) + 8*pow2(sbeta))*pow6(
        mU3)) - mQ32*(-(mA2*pow2(cbeta)*pow2(Yt)*pow4(mU3)) + pow4(Mu)*(-(mA2*
        pow2(cbeta)*pow2(Yt)) + mU32*pow2(cbeta)*pow2(Yt) - 2*(-7 + 5*pow2(
        cbeta) + 8*pow2(sbeta))*pow4(mU3)) + (mU32*(-13 + 5*pow2(cbeta) + 8*
        pow2(sbeta)) + pow2(cbeta)*pow2(Yt))*pow6(Mu) + 3*pow2(cbeta)*pow2(Yt)*
        pow6(mU3) + Mu2*(2*mA2*mU32*pow2(cbeta)*pow2(Yt) - 5*pow2(cbeta)*pow2(
        Yt)*pow4(mU3) + (-9 + 5*pow2(cbeta) + 8*pow2(sbeta))*pow6(mU3)) + 2*
        pow8(mU3))))/(mQ32*(mQ32 - Mu2)*mU32*pow2(Mu2 - mU32)*pow2(sbeta)) + 6*
        lmQ3MR*((-10*Xt*Yt*pow2(cbeta))/((mQ32 - mU32)*pow2(sbeta)) + (2*(-6*
        mA2 + 7*mQ32 + 5*mU32)*Yt*pow2(cbeta)*pow3(Xt))/(pow2(sbeta)*pow3(mQ32
        - mU32)) + (6*mQ32*Xt6)/pow4(mQ32 - mU32) - (Xt2*(pow2(cbeta)*pow2(mQ32
        - Mu2)*pow2(Mu2 - mU32)*(mA2 + 3*mQ32 - 5*mU32 + pow2(Yt)) - 2*mU32*
        pow4(Mu)*(-2*Mu2*mU32*(7 + 9*pow2(sbeta)) + (8 + 9*pow2(sbeta))*pow4(
        Mu) + (7 + 9*pow2(sbeta))*pow4(mU3)) + 2*(-4*Mu2*mU32*(3 + 4*pow2(
        sbeta)) + (7 + 8*pow2(sbeta))*pow4(Mu) + 2*(3 + 4*pow2(sbeta))*pow4(
        mU3))*pow6(mQ3) - 2*pow4(mQ3)*(-(mU32*(18 + 23*pow2(sbeta))*pow4(Mu)) -
        2*Mu2*(-1 + pow2(sbeta))*pow4(mU3) + 2*(7 + 8*pow2(sbeta))*pow6(Mu) + (
        5 + 9*pow2(sbeta))*pow6(mU3)) + 2*mQ32*Mu2*(2*mU32*(-2 + pow2(sbeta))*
        pow4(Mu) - 4*Mu2*(3 + 7*pow2(sbeta))*pow4(mU3) + (9 + 8*pow2(sbeta))*
        pow6(Mu) + 2*(5 + 9*pow2(sbeta))*pow6(mU3))))/(pow2(mQ32 - Mu2)*pow2(
        mQ32 - mU32)*pow2(Mu2 - mU32)*pow2(sbeta)) + (4*pow2(cbeta)*pow2(mQ32 -
        Mu2)*pow2(Mu2 - mU32) - 2*mQ32*Mu2*mU32*(mQ32 + mU32)*(-3 + 4*pow2(
        sbeta)) + (-3 + 4*pow2(sbeta))*pow4(mQ3)*pow4(mU3) + pow4(Mu)*(4*mQ32*
        mU32*(-3 + 4*pow2(sbeta)) + 4*pow2(sbeta)*pow4(mQ3) + (-1 + 4*pow2(
        sbeta))*pow4(mU3)) + (mU32*(2 - 8*pow2(sbeta)) - 8*mQ32*pow2(sbeta))*
        pow6(Mu) + (2 + 4*pow2(sbeta))*pow8(Mu))/(pow2(mQ32 - Mu2)*pow2(Mu2 -
        mU32)*pow2(sbeta)) + (Xt4*(pow6(mQ3)*((mU32*(70 + 44*pow2(sbeta)) + 9*
        pow2(cbeta)*pow2(Yt))*pow4(Mu) + 9*pow2(cbeta)*pow2(Yt)*pow4(mU3) - 2*
        Mu2*(9*mU32*pow2(cbeta)*pow2(Yt) + (19 + 14*pow2(sbeta))*pow4(mU3)) -
        4*(9 + 5*pow2(sbeta))*pow6(Mu) + (2 + 4*pow2(sbeta))*pow6(mU3)) + (Mu2
        - mU32)*pow4(mQ3)*(3*mA2*mU32*pow2(cbeta)*pow2(Yt) - 2*(mU32*(11 + 9*
        pow2(sbeta)) + 9*pow2(cbeta)*pow2(Yt))*pow4(Mu) + 3*pow2(cbeta)*pow2(
        Yt)*pow4(mU3) - 3*Mu2*(mA2*pow2(cbeta)*pow2(Yt) - 5*mU32*pow2(cbeta)*
        pow2(Yt) + (3 + 2*pow2(sbeta))*pow4(mU3)) + 10*(2 + pow2(sbeta))*pow6(
        Mu) + (13 + 14*pow2(sbeta))*pow6(mU3)) + (-2*Mu2*mU32*(17 + 10*pow2(
        sbeta)) + 2*(9 + 5*pow2(sbeta))*pow4(Mu) + (17 + 10*pow2(sbeta))*pow4(
        mU3))*pow8(mQ3) - pow4(Mu)*(3*mA2*pow2(cbeta)*pow2(Yt)*pow4(mU3) +
        pow4(Mu)*(3*mA2*pow2(cbeta)*pow2(Yt) + 3*mU32*pow2(cbeta)*pow2(Yt) + 2*
        (8 + 7*pow2(sbeta))*pow4(mU3)) + 3*pow2(cbeta)*pow2(Yt)*pow6(mU3) - 2*
        Mu2*(3*mA2*mU32*pow2(cbeta)*pow2(Yt) + 3*pow2(cbeta)*pow2(Yt)*pow4(mU3)
        + (15 + 14*pow2(sbeta))*pow6(mU3)) + (15 + 14*pow2(sbeta))*pow8(mU3)) +
        mQ32*Mu2*(6*mA2*pow2(cbeta)*pow2(Yt)*pow4(mU3) + 2*pow4(Mu)*(3*mA2*
        pow2(cbeta)*pow2(Yt) - 6*mU32*pow2(cbeta)*pow2(Yt) + 2*(6 + 5*pow2(
        sbeta))*pow4(mU3)) + (mU32*(2 + 4*pow2(sbeta)) + 9*pow2(cbeta)*pow2(Yt)
        )*pow6(Mu) + 6*pow2(cbeta)*pow2(Yt)*pow6(mU3) - Mu2*(12*mA2*mU32*pow2(
        cbeta)*pow2(Yt) + 3*pow2(cbeta)*pow2(Yt)*pow4(mU3) + (50 + 52*pow2(
        sbeta))*pow6(mU3)) + (26 + 28*pow2(sbeta))*pow8(mU3))))/(pow2(mQ32 -
        Mu2)*pow2(Mu2 - mU32)*pow2(sbeta)*pow4(mQ32 - mU32))) + (3*Xt4*(-((
        pow2(cbeta)*(7*mU32*pow2(Yt)*pow4(mQ3) + mA2*(2*mQ32*mU32*pow2(Yt) + (
        12*mU32 - pow2(Yt))*pow4(mQ3) - pow2(Yt)*pow4(mU3)) - 3*(mU32 - pow2(
        Yt))*pow6(mQ3) + pow2(Yt)*pow6(mU3) + mQ32*(13*pow2(Yt)*pow4(mU3) + 3*
        pow6(mU3))))/mU32) - (-2*(-1 + 2*pow2(sbeta))*pow2(-(Mu*mU32) + pow3(
        Mu))*pow6(mU3) + pow6(mQ3)*(2*mU32*(113 + 50*pow2(sbeta))*pow4(Mu) -
        Mu2*(233 + 164*pow2(sbeta))*pow4(mU3) - 3*(21 + 4*pow2(sbeta))*pow6(Mu)
        + 2*(41 + 38*pow2(sbeta))*pow6(mU3)) + (-24*Mu2*mU32*(4 + pow2(sbeta))
        + (43 + 12*pow2(sbeta))*pow4(Mu) + (49 + 12*pow2(sbeta))*pow4(mU3))*
        pow8(mQ3) + mQ32*mU32*(-2*(83 + 90*pow2(sbeta))*pow4(Mu)*pow4(mU3) +
        mU32*(67 + 92*pow2(sbeta))*pow6(Mu) + 3*Mu2*(31 + 28*pow2(sbeta))*pow6(
        mU3) + 12*pow8(Mu) + 2*(-1 + 2*pow2(sbeta))*pow8(mU3)) + pow4(mQ3)*(3*(
        39 + 20*pow2(sbeta))*pow4(Mu)*pow4(mU3) - 2*mU32*(75 + 38*pow2(sbeta))*
        pow6(Mu) + 18*Mu2*(5 + 6*pow2(sbeta))*pow6(mU3) + 24*pow8(Mu) - (93 +
        92*pow2(sbeta))*pow8(mU3)))/((mQ32 - Mu2)*pow2(Mu2 - mU32))))/(mQ32*
        pow2(sbeta)*pow4(mQ32 - mU32)) + (6*lmUMR*(-(pow4(mQ3)*(4*Mu2*mU32 +
        pow4(Mu) - 2*pow4(mU3))) + 2*mQ32*(4*mU32*pow4(Mu) - 2*Mu2*pow4(mU3) +
        pow6(Mu)) - 3*pow8(Mu) + (Xt4*(pow6(mQ3)*(44*mU32*pow4(Mu) - 28*Mu2*
        pow4(mU3) - 18*pow6(Mu) + 4*pow6(mU3)) + (-8*Mu2*mU32 + 3*pow4(Mu) + 4*
        pow4(mU3))*pow8(mQ3) + pow4(mU3)*pow8(Mu) - 2*mQ32*Mu2*(9*pow4(Mu)*
        pow4(mU3) - 14*mU32*pow6(Mu) + 2*Mu2*pow6(mU3) + 6*pow8(Mu) - 2*pow8(
        mU3)) + pow4(mQ3)*(41*pow4(Mu)*pow4(mU3) - 60*mU32*pow6(Mu) - 4*Mu2*
        pow6(mU3) + 25*pow8(Mu) - 2*pow8(mU3))))/pow4(mQ32 - mU32) - (2*Xt2*((-
        4*Mu2*mU32 + pow4(Mu) + 2*pow4(mU3))*pow6(mQ3) + pow4(mQ3)*(13*mU32*
        pow4(Mu) - 6*Mu2*pow4(mU3) - 4*pow6(Mu)) - 6*pow4(mU3)*pow6(Mu) + 2*
        pow4(Mu)*pow6(mU3) + 7*mU32*pow8(Mu) + mQ32*(4*pow4(Mu)*pow4(mU3) - 10*
        mU32*pow6(Mu) + 3*pow8(Mu)) - 2*power10(Mu)))/pow2(mQ32 - mU32)))/(
        pow2(mQ32 - Mu2)*pow2(Mu2 - mU32)*pow2(sbeta))) + (6*Xt4*(-6*mQ32*(mQ32
        - Mu2)*(Mu2 - mU32)*mU32*dilog(1 - mQ32/mU32)*pow2(mQ32 - mU32)*pow2(
        sbeta) + 3*mA2*mU32*pow2(cbeta)*pow4(mQ3)*pow4(Mu) + 9*mU32*pow2(cbeta)
        *pow2(Yt)*pow4(mQ3)*pow4(Mu) + 2*mQ32*(mQ32 - Mu2)*(Mu2 - mU32)*mU32*
        dilog(1 - mQ32/Mu2)*(mQ32*(-6*Mu2 + 2*mU32) + 2*pow4(mQ3) + 3*pow4(Mu)
        - pow4(mU3)) + 3*mA2*Mu2*pow2(cbeta)*pow4(mQ3)*pow4(mU3) + 3*Mu2*pow2(
        cbeta)*pow2(Yt)*pow4(mQ3)*pow4(mU3) - 6*mA2*mQ32*pow2(cbeta)*pow4(Mu)*
        pow4(mU3) - 12*mQ32*pow2(cbeta)*pow2(Yt)*pow4(Mu)*pow4(mU3) - 73*pow4(
        mQ3)*pow4(Mu)*pow4(mU3) + 22*dilog(1 - mU32/Mu2)*pow4(mQ3)*pow4(Mu)*
        pow4(mU3) + 6*pow2(cbeta)*pow4(mQ3)*pow4(Mu)*pow4(mU3) - 114*pow2(
        sbeta)*pow4(mQ3)*pow4(Mu)*pow4(mU3) - 5*mA2*Mu2*mU32*pow2(cbeta)*pow6(
        mQ3) - 11*Mu2*mU32*pow2(cbeta)*pow2(Yt)*pow6(mQ3) + 56*mU32*pow4(Mu)*
        pow6(mQ3) + 16*mU32*dilog(1 - mU32/Mu2)*pow4(Mu)*pow6(mQ3) + 2*mA2*
        pow2(cbeta)*pow4(Mu)*pow6(mQ3) - 3*mU32*pow2(cbeta)*pow4(Mu)*pow6(mQ3)
        + 62*mU32*pow2(sbeta)*pow4(Mu)*pow6(mQ3) + 2*pow2(cbeta)*pow2(Yt)*pow4(
        Mu)*pow6(mQ3) + 21*Mu2*pow4(mU3)*pow6(mQ3) - 20*Mu2*dilog(1 - mU32/Mu2)
        *pow4(mU3)*pow6(mQ3) + 3*mA2*pow2(cbeta)*pow4(mU3)*pow6(mQ3) - 3*Mu2*
        pow2(cbeta)*pow4(mU3)*pow6(mQ3) + 52*Mu2*pow2(sbeta)*pow4(mU3)*pow6(
        mQ3) + 9*pow2(cbeta)*pow2(Yt)*pow4(mU3)*pow6(mQ3) - 12*mU32*pow4(mQ3)*
        pow6(Mu) - 18*mU32*dilog(1 - mU32/Mu2)*pow4(mQ3)*pow6(Mu) + 18*mQ32*
        pow4(mU3)*pow6(Mu) - 6*mQ32*dilog(1 - mU32/Mu2)*pow4(mU3)*pow6(Mu) - 4*
        pow6(mQ3)*pow6(Mu) + 5*mA2*mQ32*Mu2*pow2(cbeta)*pow6(mU3) + 11*mQ32*
        Mu2*pow2(cbeta)*pow2(Yt)*pow6(mU3) + 61*Mu2*pow4(mQ3)*pow6(mU3) - 2*
        Mu2*dilog(1 - mU32/Mu2)*pow4(mQ3)*pow6(mU3) - 6*mA2*pow2(cbeta)*pow4(
        mQ3)*pow6(mU3) - 3*Mu2*pow2(cbeta)*pow4(mQ3)*pow6(mU3) + 56*Mu2*pow2(
        sbeta)*pow4(mQ3)*pow6(mU3) - 12*pow2(cbeta)*pow2(Yt)*pow4(mQ3)*pow6(
        mU3) + 8*mQ32*pow4(Mu)*pow6(mU3) - 2*mQ32*dilog(1 - mU32/Mu2)*pow4(Mu)*
        pow6(mU3) + mA2*pow2(cbeta)*pow4(Mu)*pow6(mU3) - 3*mQ32*pow2(cbeta)*
        pow4(Mu)*pow6(mU3) + 58*mQ32*pow2(sbeta)*pow4(Mu)*pow6(mU3) + pow2(
        cbeta)*pow2(Yt)*pow4(Mu)*pow6(mU3) - 79*pow6(mQ3)*pow6(mU3) + 4*dilog(1
        - mU32/Mu2)*pow6(mQ3)*pow6(mU3) + 6*pow2(cbeta)*pow6(mQ3)*pow6(mU3) -
        114*pow2(sbeta)*pow6(mQ3)*pow6(mU3) - 2*pow6(Mu)*pow6(mU3) - 50*Mu2*
        mU32*pow8(mQ3) - 4*Mu2*mU32*dilog(1 - mU32/Mu2)*pow8(mQ3) - 2*mA2*Mu2*
        pow2(cbeta)*pow8(mQ3) + 2*mA2*mU32*pow2(cbeta)*pow8(mQ3) + 3*Mu2*mU32*
        pow2(cbeta)*pow8(mQ3) - 58*Mu2*mU32*pow2(sbeta)*pow8(mQ3) - 2*Mu2*pow2(
        cbeta)*pow2(Yt)*pow8(mQ3) + 2*mU32*pow2(cbeta)*pow2(Yt)*pow8(mQ3) + 6*
        pow4(Mu)*pow8(mQ3) - 4*pow2(sbeta)*pow4(Mu)*pow8(mQ3) + 46*pow4(mU3)*
        pow8(mQ3) + 4*dilog(1 - mU32/Mu2)*pow4(mU3)*pow8(mQ3) - 3*pow2(cbeta)*
        pow4(mU3)*pow8(mQ3) + 62*pow2(sbeta)*pow4(mU3)*pow8(mQ3) + 6*mQ32*mU32*
        dilog(1 - mU32/Mu2)*pow8(Mu) - 29*mQ32*Mu2*pow8(mU3) + 2*mQ32*Mu2*
        dilog(1 - mU32/Mu2)*pow8(mU3) + mA2*mQ32*pow2(cbeta)*pow8(mU3) - mA2*
        Mu2*pow2(cbeta)*pow8(mU3) + 3*mQ32*Mu2*pow2(cbeta)*pow8(mU3) - 56*mQ32*
        Mu2*pow2(sbeta)*pow8(mU3) + mQ32*pow2(cbeta)*pow2(Yt)*pow8(mU3) - Mu2*
        pow2(cbeta)*pow2(Yt)*pow8(mU3) + 30*pow4(mQ3)*pow8(mU3) - 2*dilog(1 -
        mU32/Mu2)*pow4(mQ3)*pow8(mU3) - 3*pow2(cbeta)*pow4(mQ3)*pow8(mU3) + 58*
        pow2(sbeta)*pow4(mQ3)*pow8(mU3) + 3*pow4(Mu)*pow8(mU3) - 2*pow2(sbeta)*
        pow4(Mu)*pow8(mU3) - 2*Mu2*power10(mQ3) + 2*mU32*power10(mQ3) + 4*Mu2*
        pow2(sbeta)*power10(mQ3) - 4*mU32*pow2(sbeta)*power10(mQ3) + mQ32*
        power10(mU3) - Mu2*power10(mU3) - 2*mQ32*pow2(sbeta)*power10(mU3) + 2*
        Mu2*pow2(sbeta)*power10(mU3)))/(mQ32*(mQ32 - Mu2)*(Mu2 - mU32)*mU32*
        pow2(sbeta)*pow4(mQ32 - mU32)) + 3*lmQ3MR*((48*Yt*pow2(cbeta)*pow3(Xt))
        /(pow2(mQ32 - mU32)*pow2(sbeta)) + (pow2(cbeta)*pow2(-mQ32 + mU32 +
        Xt2)*pow2(Yt)*(pow2(mQ32 - mU32) - pow4(mA)))/(mU32*deltaxyz(mA2,mU32,
        mQ32)*pow2(mQ32 - mU32)*pow2(sbeta)) + (2*Xt6*(5*mQ32*mU32 + 2*pow4(
        mQ3) - pow4(mU3)))/(mQ32*mU32*pow3(mQ32 - mU32)) + (pow2(cbeta)*pow2(
        Yt)*(mA2*(3*mQ32 + 5*mU32)*pow2(mQ32 - mU32) - pow4(mA)*(4*mQ32*mU32 +
        3*pow4(mQ3) + 5*pow4(mU3)) - pow4(mQ32 - mU32) + (mQ32 + mU32)*pow6(mA)
        + (2*Xt2*(pow4(mA)*(2*mQ32*mU32 + 3*pow4(mQ3) - pow4(mU3)) + pow4(mQ32
        - mU32) - mQ32*pow6(mA) + mA2*(2*mU32*pow4(mQ3) + mQ32*pow4(mU3) - 3*
        pow6(mQ3))))/(mQ32 - mU32) + (Xt4*(-3*(mQ32 + mU32)*pow4(mA) + 3*mU32*
        pow4(mQ3) + mQ32*pow4(mU3) + mA2*(3*pow4(mQ3) + 5*pow4(mU3)) + pow6(mA)
        - pow6(mQ3) - 3*pow6(mU3)))/(mQ32 - mU32)))/(mQ32*deltaxyz(mA2,mQ32,
        mU32)*pow2(sbeta)*pow4(mU3)) - (8*lmUMR*(2*mQ32*Mu2*mU32*(mQ32 + mU32)
        - pow4(mQ3)*pow4(mU3) + pow4(Mu)*(-4*mQ32*mU32 + pow4(mU3)) - 2*mU32*
        pow6(Mu) + 2*pow8(Mu)))/(pow2(mQ32 - Mu2)*pow2(Mu2 - mU32)*pow2(sbeta))
        + (2*Xt2*(Mu2*(Mu2 - mU32)*(mU32*(-2 + pow2(sbeta)) - pow2(cbeta)*pow2(
        Yt))*pow4(mU3) + (Mu2 - mU32)*(2*mU32*(-2 + pow2(sbeta)) - pow2(cbeta)*
        pow2(Yt))*pow6(mQ3) + pow4(mQ3)*(-(mA2*mU32*pow2(cbeta)*pow2(Yt)) + (-
        2*mU32*(-2 + pow2(sbeta)) + pow2(cbeta)*pow2(Yt))*pow4(Mu) + pow2(
        cbeta)*pow2(Yt)*pow4(mU3) + Mu2*(mA2*pow2(cbeta)*pow2(Yt) - 2*mU32*
        pow2(cbeta)*pow2(Yt) + (-12 + 2*pow2(cbeta) + 3*pow2(sbeta))*pow4(mU3))
        - (-4 + 2*pow2(cbeta) + pow2(sbeta))*pow6(mU3)) + mQ32*(-(pow4(Mu)*(
        mA2*pow2(cbeta)*pow2(Yt) - mU32*pow2(cbeta)*pow2(Yt) + (-4 + 2*pow2(
        cbeta) + pow2(sbeta))*pow4(mU3))) - pow2(cbeta)*pow2(Yt)*pow6(mU3) +
        Mu2*(mA2*mU32*pow2(cbeta)*pow2(Yt) + 2*(1 + pow2(cbeta))*pow6(mU3)) + (
        -2 + pow2(sbeta))*pow8(mU3))))/(mQ32*(mQ32 - Mu2)*(mQ32 - mU32)*(Mu2 -
        mU32)*pow2(sbeta)*pow4(mU3)) - (Xt4*(-(pow2(Mu2 - mU32)*pow4(Mu)*pow4(
        mU3)*(-(mA2*pow2(cbeta)*pow2(Yt)) + mU32*pow2(cbeta)*pow2(Yt) + (-2 +
        4*pow2(sbeta))*pow4(mU3))) + pow8(mQ3)*(mA2*pow2(cbeta)*pow2(Yt)*pow4(
        mU3) + pow4(Mu)*(mA2*pow2(cbeta)*pow2(Yt) - 4*mU32*pow2(cbeta)*pow2(Yt)
        - 2*(47 + 18*pow2(sbeta))*pow4(mU3)) + 2*(mU32*(4 - 8*pow2(sbeta)) +
        pow2(cbeta)*pow2(Yt))*pow6(Mu) + 2*Mu2*(-(mA2*mU32*pow2(cbeta)*pow2(Yt)
        ) + pow2(cbeta)*pow2(Yt)*pow4(mU3) + (78 + 60*pow2(sbeta))*pow6(mU3)) -
        2*(37 + 34*pow2(sbeta))*pow8(mU3)) + pow6(mQ3)*(2*(-(mA2*pow2(cbeta)*
        pow2(Yt)) + mU32*pow2(cbeta)*pow2(Yt) + (88 + 60*pow2(sbeta))*pow4(mU3)
        )*pow6(Mu) + 2*(-(mA2*pow2(cbeta)*pow2(Yt)) - 5*mU32*pow2(cbeta)*pow2(
        Yt) + (26 + 32*pow2(sbeta))*pow4(mU3))*pow6(mU3) - pow4(Mu)*(-2*mA2*
        mU32*pow2(cbeta)*pow2(Yt) + 11*pow2(cbeta)*pow2(Yt)*pow4(mU3) + 4*(67 +
        50*pow2(sbeta))*pow6(mU3)) + (mU32*(-4 + 8*pow2(sbeta)) - pow2(cbeta)*
        pow2(Yt))*pow8(Mu) + 2*Mu2*(mA2*pow2(cbeta)*pow2(Yt)*pow4(mU3) + 10*
        pow2(cbeta)*pow2(Yt)*pow6(mU3) + 4*(7 + pow2(sbeta))*pow8(mU3))) +
        pow4(mQ3)*(-2*Mu2*(-(mA2*pow2(cbeta)*pow2(Yt)) - 11*mU32*pow2(cbeta)*
        pow2(Yt) + 6*(9 + 10*pow2(sbeta))*pow4(mU3))*pow6(mU3) + 2*pow6(Mu)*(
        mA2*mU32*pow2(cbeta)*pow2(Yt) + 10*pow2(cbeta)*pow2(Yt)*pow4(mU3) + (50
        + 4*pow2(sbeta))*pow6(mU3)) + (mA2*pow2(cbeta)*pow2(Yt) - 2*(55 + 34*
        pow2(sbeta))*pow4(mU3))*pow8(Mu) + (mA2*pow2(cbeta)*pow2(Yt) - mU32*
        pow2(cbeta)*pow2(Yt) + (2 - 4*pow2(sbeta))*pow4(mU3))*pow8(mU3) + pow4(
        Mu)*(-6*mA2*pow2(cbeta)*pow2(Yt)*pow4(mU3) - 41*pow2(cbeta)*pow2(Yt)*
        pow6(mU3) + 8*(13 + 23*pow2(sbeta))*pow8(mU3))) + 2*mQ32*Mu2*mU32*((-(
        mA2*pow2(cbeta)*pow2(Yt)) - 5*mU32*pow2(cbeta)*pow2(Yt) + 4*(5 + 8*
        pow2(sbeta))*pow4(mU3))*pow6(Mu) + (-(mA2*pow2(cbeta)*pow2(Yt)) + mU32*
        pow2(cbeta)*pow2(Yt) + (-2 + 4*pow2(sbeta))*pow4(mU3))*pow6(mU3) +
        pow4(Mu)*(mA2*mU32*pow2(cbeta)*pow2(Yt) + 11*pow2(cbeta)*pow2(Yt)*pow4(
        mU3) - 4*(14 + 15*pow2(sbeta))*pow6(mU3)) + 6*mU32*pow8(Mu) + Mu2*(mA2*
        pow2(cbeta)*pow2(Yt)*pow4(mU3) - 7*pow2(cbeta)*pow2(Yt)*pow6(mU3) + (34
        + 24*pow2(sbeta))*pow8(mU3))) + pow2(Mu2 - mU32)*(mU32*(-4 + 8*pow2(
        sbeta)) - pow2(cbeta)*pow2(Yt))*power10(mQ3)))/(mQ32*pow2(mQ32 - Mu2)*
        pow2(Mu2 - mU32)*pow2(sbeta)*pow3(mQ32 - mU32)*pow4(mU3)) + (mU32*pow2(
        Mu2 - mU32)*pow4(Mu)*(-(mA2*pow2(cbeta)*pow2(Yt)) + 3*mU32*pow2(cbeta)*
        pow2(Yt) + 2*pow4(mU3)) + pow2(Mu2 - mU32)*(4*mU32 + pow2(cbeta)*pow2(
        Yt))*pow8(mQ3) - pow6(mQ3)*(mA2*pow2(cbeta)*pow2(Yt)*pow4(mU3) + pow4(
        Mu)*(mA2*pow2(cbeta)*pow2(Yt) - 5*mU32*pow2(cbeta)*pow2(Yt) + 2*(-16 +
        7*pow2(cbeta) + 10*pow2(sbeta))*pow4(mU3)) + 2*(4*mU32 + pow2(cbeta)*
        pow2(Yt))*pow6(Mu) - pow2(cbeta)*pow2(Yt)*pow6(mU3) - 2*Mu2*(mA2*mU32*
        pow2(cbeta)*pow2(Yt) - 2*pow2(cbeta)*pow2(Yt)*pow4(mU3) + 2*(-6 + 7*
        pow2(cbeta) + 10*pow2(sbeta))*pow6(mU3)) + 2*(-3 + 7*pow2(cbeta) + 10*
        pow2(sbeta))*pow8(mU3)) + pow4(mQ3)*(2*(mA2*pow2(cbeta)*pow2(Yt) - 2*
        mU32*pow2(cbeta)*pow2(Yt) + 2*(-10 + 7*pow2(cbeta) + 10*pow2(sbeta))*
        pow4(mU3))*pow6(Mu) + (-(mA2*pow2(cbeta)*pow2(Yt)) + 3*mU32*pow2(cbeta)
        *pow2(Yt) + 2*pow4(mU3))*pow6(mU3) + pow4(Mu)*(-5*mA2*mU32*pow2(cbeta)*
        pow2(Yt) + 8*pow2(cbeta)*pow2(Yt)*pow4(mU3) + (38 - 56*pow2(cbeta) -
        80*pow2(sbeta))*pow6(mU3)) + (4*mU32 + pow2(cbeta)*pow2(Yt))*pow8(Mu) +
        4*Mu2*(mA2*pow2(cbeta)*pow2(Yt)*pow4(mU3) - 2*pow2(cbeta)*pow2(Yt)*
        pow6(mU3) + (-4 + 7*pow2(cbeta) + 10*pow2(sbeta))*pow8(mU3))) - mQ32*
        Mu2*((mA2*pow2(cbeta)*pow2(Yt) - mU32*pow2(cbeta)*pow2(Yt) + 2*(-14 +
        7*pow2(cbeta) + 10*pow2(sbeta))*pow4(mU3))*pow6(Mu) - 2*mA2*pow2(cbeta)
        *pow2(Yt)*pow6(mU3) - 4*pow4(Mu)*(mA2*mU32*pow2(cbeta)*pow2(Yt) - 2*
        pow2(cbeta)*pow2(Yt)*pow4(mU3) + (-11 + 7*pow2(cbeta) + 10*pow2(sbeta))
        *pow6(mU3)) + 6*pow2(cbeta)*pow2(Yt)*pow8(mU3) + Mu2*(5*mA2*pow2(cbeta)
        *pow2(Yt)*pow4(mU3) - 13*pow2(cbeta)*pow2(Yt)*pow6(mU3) + 2*(-13 + 7*
        pow2(cbeta) + 10*pow2(sbeta))*pow8(mU3)) + 4*power10(mU3)))/(mQ32*pow2(
        mQ32 - Mu2)*pow2(Mu2 - mU32)*pow2(sbeta)*pow4(mU3)));
}

double ThresholdCalculator::getDeltaLambdaYtau6(int omitLogs) const
{
   using std::log;
   using himalaya::dilog;

   const double MR2 = pow2(p.scale);
   const double lMR = omitLogs*log(MR2);
   const double mL32 = p.ml2(2,2);
   const double mE32 = p.me2(2,2);
   const double mL3 = std::sqrt(mL32);
   const double mE3 = std::sqrt(mE32);
   const double mA = p.MA;
   const double beta = std::atan(p.vu/p.vd);
   const double sbeta = std::sin(beta);
   const double cbeta = std::cos(beta);
   const double sbe = sbeta;
   const double cbe = cbeta;
   const double Xtau = p.Ae(2,2) - p.mu*p.vu/p.vd;
   const double Ytau = p.Ae(2,2) + p.mu*p.vd/p.vu;
   const double Mu = p.mu;

   if(std::abs(mE3 - mL3) < 0.1*mE3){
     const double lMS = log(mL32/MR2);
     const double M = mL3;
     const double Xtau2 = pow2(Xtau);
     return (pow2(cbe)*(cbe*(4.994969600000008 - 48.*lMS)*sbe*pow3(M)*pow3(Xtau)
     + Xtau2*pow4(M)*(-0.1268863999999965 + pow2(cbe)*(lMS*(99. - 40.5*pow2(sbe))
     + 25.88317440000001*pow2(sbe)) + (-96. + 6.75*lMS)*pow4(cbe) + lMS*(74.25
     - 99.*pow2(sbe) + 6.75*pow4(sbe))) + pow2(M)*(-0.12437120000000013
     + pow2(cbe)*(6.621856000000005*pow2(sbe) + lMS*(-22.5 + 4.5*pow2(sbe)))
     + (39. - 0.75*lMS)*pow4(cbe) + lMS*(-18.75 + 22.5*pow2(sbe)
     - 0.75*pow4(sbe)))*pow4(Xtau) + cbe*(4.507545599999993 + 108.*lMS)*sbe
     *Xtau*pow5(M) + cbe*(-0.24874240000000025 + 6.*lMS)*M*sbe*pow5(Xtau)
     + (0.7512575999999989 + pow2(cbe)*(3.9004104065361593*pow2(sbe)
     + lMS*(-30. + 58.5*pow2(sbe))) + (54. - 9.75*lMS + 18.*pow2(lMS))
     *pow4(cbe) + lMS*(-2.25 + 30.*pow2(sbe) - 9.75*pow4(sbe)))*pow6(M)
     + (lMS*(1.5 - 1.5*pow2(sbe)) + pow2(cbe)*(1.5*lMS - 0.12437120000000013
     *pow2(sbe)) - 3.*pow4(cbe))*pow6(Xtau)))/pow6(M);
   }
   return 2*(5 + lMR*(10 - 8*log(mE3)) - 8*log(mE3) - 12*(1 + lMR)*log(mL3) + 5*
        pow2(lMR) + 8*pow2(log(mE3)) + 12*pow2(log(mL3)) - (pow2(Xtau)*(-36*
        pow2(mE3)*(pow2(mE3) - pow2(mL3))*pow2(mL3)*pow2(log(mE3)) - 4*pow2(
        mL3)*pow4(mE3) - 2*lMR*pow2(mL3)*pow4(mE3) - 38*log(mL3)*pow2(mL3)*
        pow4(mE3) - 36*lMR*log(mL3)*pow2(mL3)*pow4(mE3) + 36*pow2(mL3)*pow2(
        log(mL3))*pow4(mE3) + pow2(mE3)*pow4(mL3) - lMR*pow2(mE3)*pow4(mL3) +
        34*log(mL3)*pow2(mE3)*pow4(mL3) + 32*lMR*log(mL3)*pow2(mE3)*pow4(mL3) -
        28*pow2(mE3)*pow2(log(mL3))*pow4(mL3) - 2*log(mE3)*pow2(mE3)*(-3*(7 +
        6*lMR)*pow2(mE3)*pow2(mL3) + pow4(mE3) + 16*(1 + lMR)*pow4(mL3) + 4*
        log(mL3)*pow4(mL3)) + dilog(1 - pow2(mL3)/pow2(mE3))*(-6*pow2(mL3)*
        pow4(mE3) + 6*pow2(mE3)*pow4(mL3)) + pow6(mE3) + lMR*pow6(mE3) + 2*
        pow6(mL3) + 2*lMR*pow6(mL3) - 4*log(mL3)*pow6(mL3)))/(pow2(mE3)*pow2(
        mL3)*pow2(pow2(mE3) - pow2(mL3))) + (2*pow4(Xtau)*(-28*pow2(mL3)*pow4(
        mE3) - 16*lMR*pow2(mL3)*pow4(mE3) - 12*log(mL3)*pow2(mL3)*pow4(mE3) -
        14*lMR*log(mL3)*pow2(mL3)*pow4(mE3) + 18*pow2(mL3)*pow2(log(mL3))*pow4(
        mE3) + 29*pow2(mE3)*pow4(mL3) + 17*lMR*pow2(mE3)*pow4(mL3) - 40*log(
        mL3)*pow2(mE3)*pow4(mL3) - 10*lMR*log(mL3)*pow2(mE3)*pow4(mL3) + 14*
        pow2(mE3)*pow2(log(mL3))*pow4(mL3) - 2*log(mE3)*pow2(mE3)*(-((22 + 7*
        lMR)*pow2(mE3)*pow2(mL3)) + 4*log(mL3)*pow2(mL3)*(pow2(mE3) + pow2(mL3)
        ) + pow4(mE3) - (3 + 5*lMR)*pow4(mL3)) + dilog(1 - pow2(mL3)/pow2(mE3))
        *(-3*pow2(mL3)*pow4(mE3) + 3*pow2(mE3)*pow4(mL3)) - 2*pow2(log(mE3))*(
        5*pow2(mL3)*pow4(mE3) + 3*pow2(mE3)*pow4(mL3)) + pow6(mE3) + lMR*pow6(
        mE3) - 2*pow6(mL3) - 2*lMR*pow6(mL3) + 4*log(mL3)*pow6(mL3)))/(pow2(
        mE3)*pow2(mL3)*pow3(pow2(mE3) - pow2(mL3))) - ((4*dilog(1 - pow2(mE3)/
        pow2(mL3))*pow2(mE3)*(pow2(mE3) - pow2(mL3))*pow2(mL3) - 12*pow2(mL3)*
        pow4(mE3) - 6*lMR*pow2(mL3)*pow4(mE3) + 6*log(mE3)*pow2(mL3)*pow4(mE3)
        + 6*log(mL3)*pow2(mL3)*pow4(mE3) - 40*log(mE3)*log(mL3)*pow2(mL3)*pow4(
        mE3) + 20*pow2(mL3)*pow2(log(mE3))*pow4(mE3) + 20*pow2(mL3)*pow2(log(
        mL3))*pow4(mE3) + 9*pow2(mE3)*pow4(mL3) + 3*lMR*pow2(mE3)*pow4(mL3) +
        20*log(mE3)*pow2(mE3)*pow4(mL3) + 12*lMR*log(mE3)*pow2(mE3)*pow4(mL3) -
        26*log(mL3)*pow2(mE3)*pow4(mL3) - 12*lMR*log(mL3)*pow2(mE3)*pow4(mL3) -
        48*log(mE3)*log(mL3)*pow2(mE3)*pow4(mL3) + 12*pow2(mE3)*pow2(log(mE3))*
        pow4(mL3) + 36*pow2(mE3)*pow2(log(mL3))*pow4(mL3) + dilog(1 - pow2(mL3)
        /pow2(mE3))*(-2*pow2(mL3)*pow4(mE3) + 2*pow2(mE3)*pow4(mL3)) + pow6(
        mE3) + lMR*pow6(mE3) - 2*log(mE3)*pow6(mE3) + 2*pow6(mL3) + 2*lMR*pow6(
        mL3) - 4*log(mL3)*pow6(mL3))*pow6(Xtau))/(pow2(mE3)*pow2(mL3)*pow4(
        pow2(mE3) - pow2(mL3))) + (pow2(sbeta)*(1.5 - 5*log(mE3) - 3*(3 + 2*
        lMR)*log(mL3) + 2*pow2(lMR) - (2*pow2(mA))/pow2(mE3) + lMR*(-2*log(mE3)
        + pow2(mA)*(-2/pow2(mE3) - 1/pow2(mL3))) + 2*log(mA)*(7 - 6*log(mE3) -
        6*log(mL3) + pow2(mA)*(2/pow2(mE3) + 1/pow2(mL3))) - pow2(mA)/pow2(mL3)
        + pow2(Pi) + 12*pow2(log(mA)) + 8*pow2(log(mE3)) + 12*pow2(log(mL3)) +
        (phixyz(pow2(mA),pow2(mE3),pow2(mE3))*(-deltaxyz(pow2(mA),pow2(mE3),
        pow2(mE3)) - 6*pow2(mA)*pow2(mE3) + pow4(mA)))/pow4(mE3) + (phixyz(
        pow2(mA),pow2(mL3),pow2(mL3))*(-deltaxyz(pow2(mA),pow2(mL3),pow2(mL3))
        - 7*pow2(mA)*pow2(mL3) + pow4(mA)))/pow4(mL3) + (2*Xtau*Ytau*(phixyz(
        pow2(mA),pow2(mE3),pow2(mE3))*(-deltaxyz(pow2(mA),pow2(mE3),pow2(mE3))
        - 6*pow2(mA)*pow2(mE3) + pow4(mA))*pow4(mL3) + pow4(mE3)*(-(phixyz(
        pow2(mA),pow2(mE3),pow2(mL3))*pow2(mL3)*(pow2(mA) - pow2(mE3) + pow2(
        mL3))) + phixyz(pow2(mA),pow2(mL3),pow2(mL3))*(deltaxyz(pow2(mA),pow2(
        mL3),pow2(mL3)) + 7*pow2(mA)*pow2(mL3) - pow4(mA)) + 4*log(mE3/mL3)*(-3
        - 3*lMR + log(mA) + 2*log(mE3) + 3*log(mL3))*pow4(mL3))))/((pow2(mE3) -
        pow2(mL3))*pow4(mE3)*pow4(mL3)) + (2*Ytau*pow3(Xtau)*(phixyz(pow2(mA),
        pow2(mE3),pow2(mE3))*(pow2(mA)*(pow2(mA) - 6*pow2(mE3))*(pow2(mE3) -
        pow2(mL3)) + deltaxyz(pow2(mA),pow2(mE3),pow2(mE3))*(-3*pow2(mE3) +
        pow2(mL3)))*pow4(mL3) + pow4(mE3)*(phixyz(pow2(mA),pow2(mL3),pow2(mL3))
        *(-(deltaxyz(pow2(mA),pow2(mL3),pow2(mL3))*(pow2(mE3) - 5*pow2(mL3))) +
        pow2(mA)*(pow2(mA) - 7*pow2(mL3))*(pow2(mE3) - pow2(mL3))) + phixyz(
        pow2(mA),pow2(mE3),pow2(mL3))*pow2(mL3)*(-2*deltaxyz(pow2(mA),pow2(mE3)
        ,pow2(mL3)) - (pow2(mE3) - pow2(mL3))*(pow2(mA) - pow2(mE3) + pow2(mL3)
        )) + 4*(-3*(2 + lMR)*(pow2(mE3) - pow2(mL3)) + log(mE3)*(3*(3 + lMR)*
        pow2(mE3) + log(mL3)*(2*pow2(mA) - pow2(mE3) - 3*pow2(mL3)) + 3*(1 +
        lMR)*pow2(mL3) + log(mA)*(-6*pow2(mA) - pow2(mE3) + pow2(mL3))) - log(
        mL3)*(3*(1 + lMR)*pow2(mE3) + 3*(3 + lMR)*pow2(mL3) + log(mA)*(-6*pow2(
        mA) - pow2(mE3) + pow2(mL3))) + 2*(pow2(mA) - pow2(mE3) - pow2(mL3))*
        pow2(log(mE3)) + (-4*pow2(mA) + 3*pow2(mE3) + 5*pow2(mL3))*pow2(log(
        mL3)))*pow4(mL3))))/(pow3(pow2(mE3) - pow2(mL3))*pow4(mE3)*pow4(mL3)) +
        (pow2(Ytau)*(-4/pow2(mE3) + lMR*(-4/pow2(mE3) - 2/pow2(mL3)) - 2/pow2(
        mL3) + (phixyz(pow2(mA),pow2(mE3),pow2(mL3))*(deltaxyz(pow2(mA),pow2(
        mE3),pow2(mL3)) - pow2(pow2(mA) - pow2(mE3) + pow2(mL3))))/(deltaxyz(
        pow2(mA),pow2(mE3),pow2(mL3))*pow2(mL3)) + (2*log(mL3)*(-pow2(mA) + 2*
        pow2(mE3) + pow2(mL3) - (pow2(mE3)*(-2*pow2(mA)*pow2(mE3) + pow4(mA) +
        pow4(mE3) - pow4(mL3)))/deltaxyz(pow2(mA),pow2(mE3),pow2(mL3)) + (pow3(
        pow2(mE3) - pow2(mL3)) - 3*(pow2(mE3) + pow2(mL3))*pow4(mA) - 3*pow2(
        mA)*(pow4(mE3) - pow4(mL3)) + pow6(mA))/deltaxyz(pow2(mA),pow2(mL3),
        pow2(mE3))))/pow4(mE3) + (2*log(mE3)*(-pow2(mA) + 3*pow2(mE3) - pow2(
        mL3) - (2*pow2(mE3)*pow2(mL3)*(pow2(mA) - pow2(mE3) + pow2(mL3)))/
        deltaxyz(pow2(mA),pow2(mE3),pow2(mL3)) + (pow3(-pow2(mE3) + pow2(mL3))
        - (5*pow2(mE3) + pow2(mL3))*pow4(mA) - pow2(mA)*(4*pow2(mE3)*pow2(mL3)
        - 5*pow4(mE3) + pow4(mL3)) + pow6(mA))/deltaxyz(pow2(mA),pow2(mL3),
        pow2(mE3))))/(pow2(mE3)*pow2(mL3)) + (2*log(mA)*(3*pow2(mE3)*pow2(mL3)
        + pow2(mA)*(pow2(mE3) + pow2(mL3)) - (pow2(mE3)*pow2(mL3)*(pow2(pow2(
        mE3) - pow2(mL3)) - pow4(mA)))/deltaxyz(pow2(mA),pow2(mE3),pow2(mL3)) -
        pow4(mE3) - pow4(mL3) + (-(pow2(mA)*(5*pow2(mE3) + 3*pow2(mL3))*pow2(
        pow2(mE3) - pow2(mL3))) + pow4(mA)*(4*pow2(mE3)*pow2(mL3) + 5*pow4(mE3)
        + 3*pow4(mL3)) + pow4(pow2(mE3) - pow2(mL3)) - (pow2(mE3) + pow2(mL3))*
        pow6(mA))/deltaxyz(pow2(mA),pow2(mL3),pow2(mE3))))/(pow2(mL3)*pow4(mE3)
        ) + (phixyz(pow2(mA),pow2(mL3),pow2(mE3))*(2*pow2(deltaxyz(pow2(mA),
        pow2(mL3),pow2(mE3))) - 4*pow2(mA)*(pow2(mE3) + pow2(mL3))*pow2(pow2(
        mE3) - pow2(mL3)) - 3*deltaxyz(pow2(mA),pow2(mL3),pow2(mE3))*(-2*pow2(
        mA)*(pow2(mE3) + pow2(mL3)) + pow2(pow2(mE3) - pow2(mL3)) + pow4(mA)) +
        2*pow4(mA)*(2*pow2(mE3)*pow2(mL3) + pow4(mE3) + 3*pow4(mL3)) + pow4(
        pow2(mE3) - pow2(mL3)) - 4*(pow2(mE3) + pow2(mL3))*pow6(mA) + pow8(mA))
        )/(deltaxyz(pow2(mA),pow2(mL3),pow2(mE3))*pow6(mE3))))/2. + pow4(Xtau)*
        ((-3*pow2(mE3) + 3*pow2(mL3) + pow2(mA)*(-5 + pow2(mE3)/pow2(mL3) - (2*
        pow2(mL3))/pow2(mE3)))/pow3(-pow2(mE3) + pow2(mL3)) - (3*log(mE3)*(4*
        pow2(mA)*pow2(mL3) + pow4(mE3) - pow4(mL3)))/pow4(pow2(mE3) - pow2(mL3)
        ) + (3*log(mL3)*(4*(1 + lMR)*pow2(mA)*pow2(mL3) + (1 + 2*lMR)*pow4(mE3)
        - (1 + 2*lMR)*pow4(mL3)))/pow4(pow2(mE3) - pow2(mL3)) + (2*log(mA)*(6*
        log(mE3)*pow2(mE3)*pow2(mL3)*(2*pow2(mA)*pow2(mL3) + pow4(mE3) - pow4(
        mL3)) - 6*log(mL3)*pow2(mE3)*pow2(mL3)*(2*pow2(mA)*pow2(mL3) + pow4(
        mE3) - pow4(mL3)) + (pow2(mE3) - pow2(mL3))*(-6*pow2(mL3)*pow4(mE3) +
        pow2(mA)*(-5*pow2(mE3)*pow2(mL3) + pow4(mE3) - 2*pow4(mL3)) + 6*pow2(
        mE3)*pow4(mL3))))/(pow2(mE3)*pow2(mL3)*pow4(pow2(mE3) - pow2(mL3))) + (
        lMR*(-6*log(mE3)*(2*pow2(mA)*pow2(mL3) + pow4(mE3) - pow4(mL3)) + ((-
        pow2(mE3) + pow2(mL3))*(-6*pow2(mL3)*pow4(mE3) + pow2(mA)*(-5*pow2(mE3)
        *pow2(mL3) + pow4(mE3) - 2*pow4(mL3)) + 6*pow2(mE3)*pow4(mL3)))/(pow2(
        mE3)*pow2(mL3))))/pow4(pow2(mE3) - pow2(mL3)) + pow2(Ytau)*((-11*pow2(
        mE3)*pow2(mL3) + pow4(mE3) - 2*pow4(mL3))/(pow2(mE3)*pow2(mL3)*pow3(-
        pow2(mE3) + pow2(mL3))) - (8*(-2*pow2(mA) + 3*pow2(mE3) + pow2(mL3))*
        pow2(log(mE3)))/pow4(pow2(mE3) - pow2(mL3)) + (4*(7*pow2(mA) - 3*pow2(
        mE3) - 11*pow2(mL3))*pow2(log(mL3)))/pow4(pow2(mE3) - pow2(mL3)) - (
        phixyz(pow2(mA),pow2(mE3),pow2(mL3))*(deltaxyz(pow2(mA),pow2(mE3),pow2(
        mL3))*(-pow2(mE3) + pow2(mL3))*(4*pow2(mA) - 3*pow2(mE3) + 3*pow2(mL3))
        - 6*pow2(deltaxyz(pow2(mA),pow2(mE3),pow2(mL3))) + pow2(pow2(mE3) -
        pow2(mL3))*pow2(pow2(mA) - pow2(mE3) + pow2(mL3))))/(2.*deltaxyz(pow2(
        mA),pow2(mE3),pow2(mL3))*pow2(mL3)*pow4(pow2(mE3) - pow2(mL3))) + (
        phixyz(pow2(mA),pow2(mE3),pow2(mE3))*(pow2(mA)*(pow2(mA) - 6*pow2(mE3))
        *(pow2(mE3) - pow2(mL3)) + deltaxyz(pow2(mA),pow2(mE3),pow2(mE3))*(-5*
        pow2(mE3) + pow2(mL3))))/pow4(-(mE3*pow2(mL3)) + pow3(mE3)) + (phixyz(
        pow2(mA),pow2(mL3),pow2(mL3))*(deltaxyz(pow2(mA),pow2(mL3),pow2(mL3))*(
        pow2(mE3) - 8*pow2(mL3)) + pow2(mA)*(pow2(mA) - 7*pow2(mL3))*(-pow2(
        mE3) + pow2(mL3))))/pow4(-(mL3*pow2(mE3)) + pow3(mL3)) - (lMR*(-6*pow2(
        mL3)*pow4(mE3) + 3*pow2(mE3)*pow4(mL3) + 12*log(mE3)*pow2(mE3)*pow4(
        mL3) + pow6(mE3) + 2*pow6(mL3)))/(pow2(mE3)*pow2(mL3)*pow4(pow2(mE3) -
        pow2(mL3))) + (log(mA)*(12*log(mE3)*(pow2(mA) + pow2(mE3) - pow2(mL3))
        - 12*log(mL3)*(pow2(mA) + pow2(mE3) - pow2(mL3)) - (pow2(pow2(mE3) -
        pow2(mL3))*(pow2(mE3)*pow2(mL3)*(pow2(pow2(mE3) - pow2(mL3)) - pow4(mA)
        ) + deltaxyz(pow2(mA),pow2(mE3),pow2(mL3))*(pow2(mA)*(pow2(mE3) - pow2(
        mL3)) - 3*pow2(mE3)*pow2(mL3) - 3*pow4(mE3) + pow4(mL3))))/(deltaxyz(
        pow2(mA),pow2(mE3),pow2(mL3))*pow2(mL3)*pow4(mE3)) + (pow3(-pow2(mE3) +
        pow2(mL3))*(3*(pow2(mE3) + pow2(mL3))*pow4(mA) - pow2(mL3)*pow4(mE3) -
        3*pow2(mE3)*pow4(mL3) - pow2(mA)*(5*pow4(mE3) + 3*pow4(mL3)) - pow6(mA)
        + 3*pow6(mE3) + pow6(mL3)))/(deltaxyz(pow2(mA),pow2(mL3),pow2(mE3))*
        pow2(mL3)*pow4(mE3))))/pow4(pow2(mE3) - pow2(mL3)) + (log(mL3)*(12*lMR*
        pow2(mL3) + 4*log(mE3)*(-11*pow2(mA) + 9*pow2(mE3) + 13*pow2(mL3)) + (
        pow2(pow2(mE3) - pow2(mL3))*(-((5*pow2(mE3) + 3*pow2(mL3))*pow4(mA)) +
        9*pow2(mL3)*pow4(mE3) + pow2(mE3)*pow4(mL3) + pow2(mA)*(4*pow2(mE3)*
        pow2(mL3) + 9*pow4(mE3) + 3*pow4(mL3)) + pow6(mA) - 9*pow6(mE3) - pow6(
        mL3)))/(deltaxyz(pow2(mA),pow2(mL3),pow2(mE3))*pow4(mE3)) + (-(pow2(mA)
        *pow2(pow2(mE3) - pow2(mL3))) + 17*pow2(mL3)*pow4(mE3) - (pow2(-(mE3*
        pow2(mL3)) + pow3(mE3))*(-2*pow2(mA)*pow2(mE3) + pow4(mA) + pow4(mE3) -
        pow4(mL3)))/deltaxyz(pow2(mA),pow2(mE3),pow2(mL3)) + 2*pow2(mE3)*pow4(
        mL3) + 4*pow6(mE3) + pow6(mL3))/pow4(mE3)))/pow4(pow2(mE3) - pow2(mL3))
        + (log(mE3)*(pow2(mA)*pow2(pow2(mE3) - pow2(mL3)) - (2*pow2(mE3)*pow2(
        mL3)*(pow2(mA) - pow2(mE3) + pow2(mL3))*pow2(pow2(mE3) - pow2(mL3)))/
        deltaxyz(pow2(mA),pow2(mE3),pow2(mL3)) - 13*pow2(mL3)*pow4(mE3) - 7*
        pow2(mE3)*pow4(mL3) - pow6(mE3) - 3*pow6(mL3) + (pow2(pow2(mE3) - pow2(
        mL3))*((3*pow2(mE3) + 5*pow2(mL3))*pow4(mA) + 5*pow2(mL3)*pow4(mE3) -
        11*pow2(mE3)*pow4(mL3) - pow2(mA)*(4*pow2(mE3)*pow2(mL3) + 5*pow4(mE3)
        + 7*pow4(mL3)) - pow6(mA) + 3*pow6(mE3) + 3*pow6(mL3)))/deltaxyz(pow2(
        mA),pow2(mL3),pow2(mE3))))/(pow2(mE3)*pow2(mL3)*pow4(pow2(mE3) - pow2(
        mL3))) + (phixyz(pow2(mA),pow2(mL3),pow2(mE3))*(2*pow2(deltaxyz(pow2(
        mA),pow2(mL3),pow2(mE3)))*(-5*pow2(mE3)*pow2(mL3) + 12*pow4(mE3) +
        pow4(mL3)) - deltaxyz(pow2(mA),pow2(mL3),pow2(mE3))*(-pow2(mE3) + pow2(
        mL3))*(3*(-3*pow2(mE3) + pow2(mL3))*pow4(mA) + 29*pow2(mL3)*pow4(mE3) +
        6*pow2(mA)*(2*pow2(mE3)*pow2(mL3) + 3*pow4(mE3) - pow4(mL3)) - 15*pow2(
        mE3)*pow4(mL3) - 17*pow6(mE3) + 3*pow6(mL3)) + pow2(pow2(mE3) - pow2(
        mL3))*(-4*pow2(mA)*(pow2(mE3) + pow2(mL3))*pow2(pow2(mE3) - pow2(mL3))
        + pow2(pow2(mE3) - pow2(mL3))*(-2*pow2(mE3)*pow2(mL3) - 3*pow4(mE3) +
        pow4(mL3)) + pow4(mA)*(4*pow2(mE3)*pow2(mL3) + 6*pow4(mE3) + 6*pow4(
        mL3)) - 4*(pow2(mE3) + pow2(mL3))*pow6(mA) + pow8(mA))))/(2.*deltaxyz(
        pow2(mA),pow2(mL3),pow2(mE3))*pow4(pow2(mE3) - pow2(mL3))*pow6(mE3))))
        + pow2(Xtau)*((2*pow2(mA)*(pow2(mE3) - 2*pow2(mL3)) + 4*pow2(mE3)*pow2(
        mL3))/(pow2(mE3)*(pow2(mE3) - pow2(mL3))*pow2(mL3)) + (8*pow2(log(mE3))
        )/(pow2(mE3) - pow2(mL3)) + (2*log(mE3)*(2*pow2(mA) - 5*pow2(mE3) +
        pow2(mL3)))/pow2(pow2(mE3) - pow2(mL3)) + (2*log(mL3)*(-2*(1 + lMR)*
        pow2(mA) + 3*pow2(mE3) + pow2(mL3) + 2*lMR*pow2(mL3) + 2*log(mE3)*(-
        pow2(mA) + pow2(mE3) + pow2(mL3))))/pow2(pow2(mE3) - pow2(mL3)) + (2*
        lMR*(pow2(mE3) + 2*log(mE3)*(pow2(mA) - pow2(mL3)) - pow2(mL3) + pow2(
        mA)*(-3 + pow2(mE3)/pow2(mL3) + (2*pow2(mL3))/pow2(mE3))))/pow2(pow2(
        mE3) - pow2(mL3)) + (4*(pow2(mA) - 3*pow2(mE3) + pow2(mL3))*pow2(log(
        mL3)))/pow2(pow2(mE3) - pow2(mL3)) + (deltaxyz(pow2(mA),pow2(mE3),pow2(
        mL3))*phixyz(pow2(mA),pow2(mE3),pow2(mL3)))/pow2(-(mL3*pow2(mE3)) +
        pow3(mL3)) + (phixyz(pow2(mA),pow2(mL3),pow2(mL3))*(deltaxyz(pow2(mA),
        pow2(mL3),pow2(mL3))*(pow2(mE3) - 2*pow2(mL3)) + pow2(mA)*(pow2(mA) -
        7*pow2(mL3))*(-pow2(mE3) + pow2(mL3))))/(pow2(pow2(mE3) - pow2(mL3))*
        pow4(mL3)) - (4*log(mA)*(log(mE3)*pow2(mE3)*(pow2(mA) + 5*pow2(mE3) -
        5*pow2(mL3))*pow2(mL3) - log(mL3)*pow2(mE3)*(pow2(mA) + 5*pow2(mE3) -
        5*pow2(mL3))*pow2(mL3) + pow2(mA)*(-3*pow2(mE3)*pow2(mL3) + pow4(mE3) +
        2*pow4(mL3))))/(pow2(mE3)*pow2(mL3)*pow2(pow2(mE3) - pow2(mL3))) + (
        phixyz(pow2(mA),pow2(mE3),pow2(mE3))*(-deltaxyz(pow2(mA),pow2(mE3),
        pow2(mE3)) - 6*pow2(mA)*pow2(mE3) + pow4(mA)))/(-(pow2(mL3)*pow4(mE3))
        + pow6(mE3)) + pow2(Ytau)*((8*pow2(log(mE3)))/pow2(pow2(mE3) - pow2(
        mL3)) + (12*pow2(log(mL3)))/pow2(pow2(mE3) - pow2(mL3)) + (phixyz(pow2(
        mA),pow2(mE3),pow2(mL3))*(deltaxyz(pow2(mA),pow2(mE3),pow2(mL3))*pow2(
        mA) - (pow2(mE3) - pow2(mL3))*pow2(pow2(mA) - pow2(mE3) + pow2(mL3))))/
        (deltaxyz(pow2(mA),pow2(mE3),pow2(mL3))*pow2(-(mL3*pow2(mE3)) + pow3(
        mL3))) + (phixyz(pow2(mA),pow2(mE3),pow2(mE3))*(-deltaxyz(pow2(mA),
        pow2(mE3),pow2(mE3)) - 6*pow2(mA)*pow2(mE3) + pow4(mA)))/(pow2(pow2(
        mE3) - pow2(mL3))*pow4(mE3)) + (phixyz(pow2(mA),pow2(mL3),pow2(mL3))*(-
        deltaxyz(pow2(mA),pow2(mL3),pow2(mL3)) - 7*pow2(mA)*pow2(mL3) + pow4(
        mA)))/(pow2(pow2(mE3) - pow2(mL3))*pow4(mL3)) + (2*lMR*(-3*pow2(mE3)*
        pow2(mL3) + 2*log(mE3)*pow2(mE3)*pow2(mL3) + pow4(mE3) + 2*pow4(mL3)))/
        (pow2(mE3)*pow2(mL3)*pow2(pow2(mE3) - pow2(mL3))) + (2*pow2(mE3) - 4*
        pow2(mL3))/(pow2(mL3)*pow4(mE3) - pow2(mE3)*pow4(mL3)) + (2*log(mE3)*(
        pow2(mE3)*pow2(mL3) + (2*pow2(mE3)*pow2(mL3)*(-pow2(mE3) + pow2(mL3))*(
        pow2(mA) - pow2(mE3) + pow2(mL3)))/deltaxyz(pow2(mA),pow2(mE3),pow2(
        mL3)) - pow4(mE3) + 2*pow4(mL3) + ((pow2(mE3) - pow2(mL3))*(-((pow2(
        mE3) - 2*pow2(mL3))*pow2(pow2(mE3) - pow2(mL3))) + (pow2(mE3) + 2*pow2(
        mL3))*pow4(mA) - 4*pow2(mA)*(2*pow2(mE3)*pow2(mL3) + pow4(mL3))))/
        deltaxyz(pow2(mA),pow2(mL3),pow2(mE3))))/(pow2(mE3)*pow2(mL3)*pow2(
        pow2(mE3) - pow2(mL3))) + (2*log(mL3)*(-2*lMR - 10*log(mE3) + (-2*pow2(
        mE3)*pow2(mL3) + pow2(mA)*(-pow2(mE3) + pow2(mL3)) + pow4(mE3) - (pow2(
        mE3)*(pow2(mE3) - pow2(mL3))*(-2*pow2(mA)*pow2(mE3) + pow4(mA) + pow4(
        mE3) - pow4(mL3)))/deltaxyz(pow2(mA),pow2(mE3),pow2(mL3)) - pow4(mL3))/
        pow4(mE3) - ((pow2(mE3) - pow2(mL3))*(pow2(-(mL3*pow2(mE3)) + pow3(mL3)
        ) + (4*pow2(mE3) + 3*pow2(mL3))*pow4(mA) - pow2(mA)*(2*pow2(mE3)*pow2(
        mL3) + 7*pow4(mE3) + 3*pow4(mL3)) - pow6(mA)))/(deltaxyz(pow2(mA),pow2(
        mL3),pow2(mE3))*pow4(mE3))))/pow2(pow2(mE3) - pow2(mL3)) + (2*log(mA)*(
        log(pow2(mL3)/pow2(mE3)) + ((-pow2(mE3) + pow2(mL3))*(deltaxyz(pow2(mA)
        ,pow2(mL3),pow2(mE3))*(pow2(mE3)*pow2(mL3)*(pow2(pow2(mE3) - pow2(mL3))
        - pow4(mA)) + deltaxyz(pow2(mA),pow2(mE3),pow2(mL3))*(-(pow2(mA)*pow2(
        mL3)) - 3*pow2(mE3)*pow2(mL3) + pow4(mE3) + pow4(mL3))) + deltaxyz(
        pow2(mA),pow2(mE3),pow2(mL3))*(pow4(mA)*(-2*pow2(mE3)*pow2(mL3) + pow4(
        mE3) - 3*pow4(mL3)) - pow2(mA)*pow2(mL3)*(2*pow2(mE3)*pow2(mL3) + pow4(
        mE3) - 3*pow4(mL3)) - pow4(pow2(mE3) - pow2(mL3)) + pow2(mL3)*pow6(mA))
        ))/(deltaxyz(pow2(mA),pow2(mE3),pow2(mL3))*deltaxyz(pow2(mA),pow2(mL3),
        pow2(mE3))*pow2(mL3)*pow4(mE3))))/pow2(pow2(mE3) - pow2(mL3)) + (
        phixyz(pow2(mA),pow2(mL3),pow2(mE3))*(-2*(-2*pow2(mE3) + pow2(mL3))*
        pow2(deltaxyz(pow2(mA),pow2(mL3),pow2(mE3))) + deltaxyz(pow2(mA),pow2(
        mL3),pow2(mE3))*((-5*pow2(mE3) + 3*pow2(mL3))*pow2(pow2(mE3) - pow2(
        mL3)) + (-5*pow2(mE3) + 3*pow2(mL3))*pow4(mA) + 2*pow2(mA)*(2*pow2(mE3)
        *pow2(mL3) + 7*pow4(mE3) - 3*pow4(mL3))) + (pow2(mE3) - pow2(mL3))*(4*
        pow2(mA)*(pow2(mE3) - pow2(mL3))*pow4(mL3) + pow4(mA)*(4*pow2(mE3)*
        pow2(mL3) + 6*pow4(mE3) + 6*pow4(mL3)) + pow4(pow2(mE3) - pow2(mL3)) -
        4*(pow2(mE3) + pow2(mL3))*pow6(mA) + pow8(mA))))/(deltaxyz(pow2(mA),
        pow2(mL3),pow2(mE3))*pow2(pow2(mE3) - pow2(mL3))*pow6(mE3))))))/pow2(
        cbeta) + (-2*pow2(lMR) + lMR*(2*log(mE3) - ((pow2(mE3) + 2*pow2(mL3))*(
        pow2(mE3) + pow2(mL3) - 2*pow2(Mu)))/(pow2(mE3)*pow2(mL3))) - (4*pow2(
        mE3)*(pow2(mE3) - 2*pow2(Mu))*pow2(log(mE3)))/pow2(pow2(mE3) - pow2(Mu)
        ) - (8*pow2(mL3)*(pow2(mL3) - 2*pow2(Mu))*pow2(log(mL3)))/pow2(pow2(
        mL3) - pow2(Mu)) + (4*dilog(1 - pow2(mL3)/pow2(Mu))*pow4(Mu))/pow2(
        pow2(mL3) - pow2(Mu)) + (2*dilog(1 - pow2(mE3)/pow2(Mu))*(2*pow2(mE3)*
        pow2(Mu) - pow4(mE3) + pow4(Mu)))/pow2(pow2(mE3) - pow2(Mu)) + (2*
        dilog(1 - pow2(mL3)/pow2(Mu))*(2*pow2(mL3)*pow2(Mu) - pow4(mL3) + pow4(
        Mu)))/pow2(pow2(mL3) - pow2(Mu)) + (-(pow2(mE3)*pow2(mL3)*(3*pow2(mE3)*
        pow2(mL3) + 2*pow4(mE3) + 4*pow4(mL3))) - 3*(7*pow2(mE3)*pow2(mL3) + 2*
        pow4(mE3) + 4*pow4(mL3))*pow4(Mu) + pow2(Mu)*(13*pow2(mL3)*pow4(mE3) +
        17*pow2(mE3)*pow4(mL3) + 2*pow6(mE3) + 4*pow6(mL3)) + 4*(pow2(mE3) + 2*
        pow2(mL3))*pow6(Mu))/(2.*pow2(mE3)*pow2(mL3)*(-pow2(mE3) + pow2(Mu))*(-
        pow2(mL3) + pow2(Mu))) + (log(mE3)*(9*pow4(mL3)*pow4(Mu) + pow4(mE3)*(-
        9*pow2(mL3)*pow2(Mu) + pow4(mL3) + 4*pow4(Mu)) + 2*(pow2(mL3) - pow2(
        Mu))*pow6(mE3) - 13*pow2(mL3)*pow6(Mu) - 2*pow2(mE3)*(2*pow2(Mu)*pow4(
        mL3) - 7*pow2(mL3)*pow4(Mu) + pow6(Mu))))/(pow2(mL3)*(pow2(mL3) - pow2(
        Mu))*pow2(pow2(mE3) - pow2(Mu))) + (log(pow2(Mu))*(2*log(mE3)*pow2(mE3)
        *pow2(mL3)*(2*pow4(mE3)*(-2*pow2(mL3)*pow2(Mu) + pow4(mL3)) - pow4(Mu)*
        (-2*pow2(mL3)*pow2(Mu) + pow4(mL3) + 3*pow4(Mu)) + pow2(mE3)*(-4*pow2(
        Mu)*pow4(mL3) + 8*pow2(mL3)*pow4(Mu))) + pow2(Mu)*(-4*pow2(-(mL3*pow2(
        Mu)) + pow3(mL3))*pow4(Mu) - 2*(pow2(mL3)*pow2(Mu) + pow4(mL3) + pow4(
        Mu))*pow6(mE3) + pow4(mE3)*(8*pow2(Mu)*pow4(mL3) + 2*pow2(mL3)*pow4(Mu)
        - 2*pow6(mL3) + 4*pow6(Mu)) + pow2(mE3)*(-8*pow4(mL3)*pow4(Mu) + 3*
        pow2(Mu)*pow6(mL3) + pow2(mL3)*pow6(Mu) - 2*pow8(Mu)))))/(pow2(mE3)*
        pow2(mL3)*pow2(pow2(mE3) - pow2(Mu))*pow2(pow2(mL3) - pow2(Mu))) + (2*
        pow2(log(pow2(Mu)))*(-2*pow2(mE3)*pow2(Mu)*(2*pow2(mL3)*pow2(Mu) -
        pow4(mL3) + pow4(Mu)) + pow4(mE3)*(2*pow2(mL3)*pow2(Mu) - pow4(mL3) +
        pow4(Mu)) + 2*pow8(Mu)))/(pow2(pow2(mE3) - pow2(Mu))*pow2(pow2(mL3) -
        pow2(Mu))) + (log(mL3)*(2*log(pow2(Mu))*pow2(mE3)*(2*pow4(mE3)*(-2*
        pow2(mL3)*pow2(Mu) + pow4(mL3) - 2*pow4(Mu)) + (-2*pow2(mL3)*pow2(Mu) +
        pow4(mL3) - 5*pow4(Mu))*pow4(Mu) + pow2(mE3)*(-4*pow2(Mu)*pow4(mL3) +
        8*pow2(mL3)*pow4(Mu) + 8*pow6(Mu))) + (pow2(mE3) - pow2(Mu))*(-4*pow2(
        mL3)*pow2(Mu)*pow2(pow2(mL3) - pow2(Mu)) + pow4(mE3)*(-6*(1 + 2*lMR)*
        pow2(mL3)*pow2(Mu) + (5 + 6*lMR)*pow4(mL3) + (13 + 6*lMR)*pow4(Mu)) +
        pow2(mE3)*(-3*(5 + 2*lMR)*pow2(Mu)*pow4(mL3) + 2*(7 + 6*lMR)*pow2(mL3)*
        pow4(Mu) + 4*pow6(mL3) - 3*(5 + 2*lMR)*pow6(Mu))) + log(mE3)*(-8*pow2(
        Mu)*pow4(mE3)*(2*pow2(mL3)*pow2(Mu) - pow4(mL3) + pow4(Mu)) + (8*pow2(
        mL3)*pow2(Mu) - 4*pow4(mL3) + 4*pow4(Mu))*pow6(mE3) + 8*pow2(mE3)*pow8(
        Mu))))/(pow2(mE3)*pow2(pow2(mE3) - pow2(Mu))*pow2(pow2(mL3) - pow2(Mu))
        ) + pow2(Xtau)*((4*dilog(1 - pow2(mE3)/pow2(Mu))*(pow2(mL3) - pow2(Mu))
        )/pow2(pow2(mE3) - pow2(mL3)) - (4*dilog(1 - pow2(mL3)/pow2(Mu))*(pow2(
        mL3) - pow2(Mu)))/pow2(pow2(mE3) - pow2(mL3)) + (2*lMR*(-2*log(mE3)*(
        pow2(mE3) - 4*pow2(mL3) + 2*pow2(Mu)) + ((pow2(mE3) - pow2(mL3))*(4*
        pow2(mL3)*pow2(Mu) - 2*pow2(mE3)*(pow2(mL3) + pow2(Mu)) + pow4(mE3) -
        2*pow4(mL3)))/(pow2(mE3)*pow2(mL3))))/pow2(pow2(mE3) - pow2(mL3)) + (4*
        pow2(mL3)*(pow2(mL3) - 2*pow2(Mu)) + pow2(mE3)*(6*pow2(mL3) + 4*pow2(
        Mu)) - 2*pow4(mE3))/(-(pow2(mL3)*pow4(mE3)) + pow2(mE3)*pow4(mL3)) + (
        8*pow2(mE3)*pow2(log(mE3))*(pow2(mE3)*(pow2(mL3) - 6*pow2(Mu)) - 2*
        pow2(mL3)*pow2(Mu) + 3*pow4(mE3) + 4*pow4(Mu)))/(pow2(pow2(mE3) - pow2(
        mL3))*pow2(pow2(mE3) - pow2(Mu))) + (8*pow2(mL3)*pow2(log(mL3))*(2*
        pow2(mE3)*(pow2(mL3) - 2*pow2(Mu)) - 6*pow2(mL3)*pow2(Mu) + 3*pow4(mL3)
        + 5*pow4(Mu)))/(pow2(pow2(mE3) - pow2(mL3))*pow2(pow2(mL3) - pow2(Mu)))
        - (2*log(mE3)*(-5*pow4(mL3)*pow4(Mu) + pow4(mE3)*(pow2(mL3)*pow2(Mu) +
        pow4(mL3) + 2*pow4(Mu)) + 2*(pow2(mL3) - pow2(Mu))*pow6(mE3) + pow2(
        mE3)*(16*pow2(Mu)*pow4(mL3) - 13*pow2(mL3)*pow4(Mu) - 11*pow6(mL3)) +
        5*pow2(Mu)*pow6(mL3) + 4*pow2(mL3)*pow6(Mu)))/((pow2(mE3) - pow2(Mu))*(
        pow2(mL3) - pow2(Mu))*pow2(-(mL3*pow2(mE3)) + pow3(mL3))) + (2*log(mL3)
        *((pow2(mE3) - pow2(Mu))*(pow2(mL3) - pow2(Mu))*(4*(pow2(mL3) - pow2(
        Mu))*pow2(Mu)*pow4(mL3) - pow4(mE3)*(-2*(1 + 5*lMR)*pow2(mL3)*pow2(Mu)
        + (11 + 8*lMR)*pow4(mL3) + (7 + 2*lMR)*pow4(Mu)) + ((7 + 2*lMR)*pow2(
        mL3) + (1 - 2*lMR)*pow2(Mu))*pow6(mE3) + pow2(mE3)*((13 + 8*lMR)*pow2(
        Mu)*pow4(mL3) - (5 + 12*lMR)*pow2(mL3)*pow4(Mu) - 4*pow6(mL3) + 4*(1 +
        lMR)*pow6(Mu))) - 4*log(mE3)*pow2(mE3)*((-10*pow2(mL3)*pow2(Mu) + 5*
        pow4(mL3) + 3*pow4(Mu))*pow6(mE3) + 3*pow4(Mu)*pow6(mL3) + 2*pow4(mE3)*
        (-9*pow2(Mu)*pow4(mL3) + 13*pow2(mL3)*pow4(Mu) + 2*pow6(mL3) - 3*pow6(
        Mu)) - 6*pow4(mL3)*pow6(Mu) + 5*pow2(mL3)*pow8(Mu) + pow2(mE3)*(22*
        pow4(mL3)*pow4(Mu) - 8*pow2(Mu)*pow6(mL3) - 24*pow2(mL3)*pow6(Mu) + 4*
        pow8(Mu))) + 2*log(pow2(Mu))*pow2(mE3)*(2*pow4(Mu)*pow6(mE3) + 2*pow4(
        mE3)*(-3*pow2(Mu)*pow4(mL3) + 2*pow2(mL3)*pow4(Mu) + pow6(mL3) - 3*
        pow6(Mu)) + pow4(Mu)*(-4*pow2(Mu)*pow4(mL3) + 3*pow2(mL3)*pow4(Mu) +
        pow6(mL3) - 2*pow6(Mu)) + pow2(mE3)*(13*pow4(mL3)*pow4(Mu) - 4*pow2(Mu)
        *pow6(mL3) - 10*pow2(mL3)*pow6(Mu) + 7*pow8(Mu)))))/(pow2(mE3)*pow2(
        pow2(mE3) - pow2(mL3))*pow2(pow2(mE3) - pow2(Mu))*pow2(pow2(mL3) -
        pow2(Mu))) - (4*log(pow2(Mu))*((pow2(mE3) - pow2(mL3))*(pow2(mE3) -
        pow2(Mu))*(pow2(mL3) - pow2(Mu))*pow2(Mu)*(pow2(Mu)*pow4(mE3) - 2*pow2(
        Mu)*pow4(mL3) + pow2(mE3)*(-(pow2(mL3)*pow2(Mu)) + pow4(mL3) - pow4(Mu)
        ) + 2*pow2(mL3)*pow4(Mu)) + log(mE3)*pow2(mE3)*pow2(mL3)*(2*pow4(Mu)*
        pow6(mE3) + 2*pow4(mE3)*(-3*pow2(Mu)*pow4(mL3) + 2*pow2(mL3)*pow4(Mu) +
        pow6(mL3) - 3*pow6(Mu)) + pow4(Mu)*(-4*pow2(Mu)*pow4(mL3) + 3*pow2(mL3)
        *pow4(Mu) + pow6(mL3) - 2*pow6(Mu)) + pow2(mE3)*(13*pow4(mL3)*pow4(Mu)
        - 4*pow2(Mu)*pow6(mL3) - 10*pow2(mL3)*pow6(Mu) + 7*pow8(Mu)))))/(pow2(
        mE3)*pow2(mL3)*pow2(pow2(mE3) - pow2(mL3))*pow2(pow2(mE3) - pow2(Mu))*
        pow2(pow2(mL3) - pow2(Mu)))) + (pow4(Xtau)*(2*dilog(1 - pow2(mE3)/pow2(
        Mu))*(-2*pow2(mE3)*pow2(mL3) + 6*pow2(mL3)*pow2(Mu) + pow4(mE3) - 2*
        pow4(mL3) - 3*pow4(Mu)) + dilog(1 - pow2(mL3)/pow2(Mu))*(4*pow2(mL3)*(
        pow2(mE3) - 3*pow2(Mu)) - 2*pow4(mE3) + 4*pow4(mL3) + 6*pow4(Mu)) - (4*
        pow2(mE3)*(pow2(mE3) + pow2(mL3))*pow2(log(mE3))*(pow2(mE3)*(pow2(mL3)
        - 14*pow2(Mu)) - 2*pow2(mL3)*pow2(Mu) + 7*pow4(mE3) + 8*pow4(Mu)))/
        pow2(pow2(mE3) - pow2(Mu)) + lMR*(2*log(mE3)*(-6*pow2(mE3)*pow2(mL3) +
        12*pow2(mL3)*pow2(Mu) + 7*pow4(mE3) - 13*pow4(mL3)) + ((-pow2(mE3) +
        pow2(mL3))*(10*pow2(mE3)*pow2(mL3)*(pow2(mE3) + pow2(Mu)) - 2*pow2(Mu)*
        pow4(mE3) + (-21*pow2(mE3) + 4*pow2(Mu))*pow4(mL3) + pow6(mE3) - 2*
        pow6(mL3)))/(pow2(mE3)*pow2(mL3))) - (8*pow2(mL3)*pow2(log(mL3))*(5*
        pow2(mE3)*pow2(pow2(mL3) - pow2(Mu)) + (pow2(mL3) - 2*pow2(Mu))*pow4(
        mE3) - 10*pow2(Mu)*pow4(mL3) + 6*pow2(mL3)*pow4(Mu) + 5*pow6(mL3)))/
        pow2(pow2(mL3) - pow2(Mu)) + ((pow2(mE3) - pow2(mL3))*(-2*pow2(Mu)*
        pow4(mL3)*(-3*pow2(mL3)*pow2(Mu) + pow4(mL3) + 2*pow4(Mu)) + (14*pow2(
        mL3)*pow2(Mu) - 15*pow4(mL3) - 3*pow4(Mu))*pow6(mE3) + pow4(mE3)*(-31*
        pow2(Mu)*pow4(mL3) + 5*pow2(mL3)*pow4(Mu) + 32*pow6(mL3) + 2*pow6(Mu))
        + (-pow2(mL3) + pow2(Mu))*pow8(mE3) + 2*pow2(mE3)*(23*pow4(mL3)*pow4(
        Mu) - 18*pow2(Mu)*pow6(mL3) - 8*pow2(mL3)*pow6(Mu) + pow8(mL3))))/(
        pow2(mE3)*pow2(mL3)*(pow2(mE3) - pow2(Mu))*(pow2(mL3) - pow2(Mu))) + (
        log(mE3)*(-2*pow2(Mu)*pow2(pow2(mE3) - pow2(Mu))*pow6(mE3) + pow6(mL3)*
        (153*pow2(Mu)*pow4(mE3) - 162*pow2(mE3)*pow4(Mu) - 50*pow6(mE3) + 47*
        pow6(Mu)) + (64*pow2(mE3)*pow2(Mu) - 33*pow4(mE3) - 27*pow4(Mu))*pow8(
        mL3) + pow4(mL3)*(-101*pow4(mE3)*pow4(Mu) - 26*pow2(Mu)*pow6(mE3) +
        118*pow2(mE3)*pow6(Mu) + 45*pow8(mE3) - 24*pow8(Mu)) + pow2(mE3)*pow2(
        mL3)*(70*pow4(mE3)*pow4(Mu) - 45*pow2(Mu)*pow6(mE3) - 19*pow2(mE3)*
        pow6(Mu) + 2*pow8(mE3) - 12*pow8(Mu))))/(pow2(mL3)*(pow2(mL3) - pow2(
        Mu))*pow2(pow2(mE3) - pow2(Mu))) - (log(mL3)*((pow2(mE3) - pow2(Mu))*(
        4*pow2(Mu)*pow2(pow2(mL3) - pow2(Mu))*pow6(mL3) + pow2(mL3)*pow4(mE3)*(
        2*(81 + 44*lMR)*pow2(Mu)*pow4(mL3) - (123 + 98*lMR)*pow2(mL3)*pow4(Mu)
        - (71 + 26*lMR)*pow6(mL3) + 4*(2 + 9*lMR)*pow6(Mu)) + pow6(mE3)*((-49 +
        10*lMR)*pow2(Mu)*pow4(mL3) + 2*(37 + 8*lMR)*pow2(mL3)*pow4(Mu) - 12*(-1
        + lMR)*pow6(mL3) - (13 + 14*lMR)*pow6(Mu)) + (-2*(25 + 14*lMR)*pow2(
        mL3)*pow2(Mu) + (27 + 14*lMR)*pow4(mL3) + (15 + 14*lMR)*pow4(Mu))*pow8(
        mE3) + pow2(mE3)*pow2(mL3)*(-2*(87 + 38*lMR)*pow4(mL3)*pow4(Mu) + (77 +
        26*lMR)*pow2(Mu)*pow6(mL3) + (145 + 74*lMR)*pow2(mL3)*pow6(Mu) - 4*
        pow8(mL3) - 12*(3 + 2*lMR)*pow8(Mu))) - 4*log(mE3)*pow2(mE3)*(2*pow4(
        mL3)*pow4(Mu)*(-10*pow2(mL3)*pow2(Mu) + 5*pow4(mL3) + 6*pow4(Mu)) + 2*
        pow6(mE3)*(-27*pow2(Mu)*pow4(mL3) + 27*pow2(mL3)*pow4(Mu) + 9*pow6(mL3)
        - 7*pow6(Mu)) + (-18*pow2(mL3)*pow2(Mu) + 9*pow4(mL3) + 7*pow4(Mu))*
        pow8(mE3) + pow4(mE3)*(95*pow4(mL3)*pow4(Mu) - 58*pow2(Mu)*pow6(mL3) -
        56*pow2(mL3)*pow6(Mu) + 11*pow8(mL3) + 8*pow8(Mu)) + pow2(mE3)*(62*
        pow4(Mu)*pow6(mL3) - 62*pow4(mL3)*pow6(Mu) - 22*pow2(Mu)*pow8(mL3) +
        18*pow2(mL3)*pow8(Mu))) - 2*log(pow2(Mu))*pow2(mE3)*(4*pow2(mL3)*(pow2(
        mL3)*pow2(Mu) - pow4(mL3) + pow4(Mu))*pow6(mE3) + 18*pow6(mL3)*pow6(Mu)
        + 2*(-2*pow2(mL3)*pow2(Mu) + pow4(mL3))*pow8(mE3) - 3*pow4(Mu)*pow8(
        mL3) - 25*pow4(mL3)*pow8(Mu) - pow4(mE3)*(41*pow4(mL3)*pow4(Mu) - 28*
        pow2(Mu)*pow6(mL3) - 18*pow2(mL3)*pow6(Mu) + 4*pow8(mL3) + pow8(Mu)) +
        4*pow2(mE3)*(-11*pow4(Mu)*pow6(mL3) + 15*pow4(mL3)*pow6(Mu) + 2*pow2(
        Mu)*pow8(mL3) - 7*pow2(mL3)*pow8(Mu)) + 12*pow2(mL3)*power10(Mu))))/(
        pow2(mE3)*pow2(pow2(mE3) - pow2(Mu))*pow2(pow2(mL3) - pow2(Mu))) - (2*
        log(pow2(Mu))*((pow2(mE3) - pow2(mL3))*pow2(Mu)*(-2*pow2(pow2(mL3) -
        pow2(Mu))*pow4(mL3)*pow4(Mu) - 2*pow6(mE3)*(pow2(mL3)*pow4(Mu) + pow6(
        mL3) + pow6(Mu)) + (pow4(mL3) + pow4(Mu))*pow8(mE3) + pow4(mE3)*(-7*
        pow4(mL3)*pow4(Mu) + 9*pow2(Mu)*pow6(mL3) + 5*pow2(mL3)*pow6(Mu) - 2*
        pow8(mL3) + pow8(Mu)) + pow2(mE3)*(-8*pow4(Mu)*pow6(mL3) + 5*pow4(mL3)*
        pow6(Mu) + 3*pow2(Mu)*pow8(mL3) - 2*pow2(mL3)*pow8(Mu))) + log(mE3)*
        pow2(mE3)*pow2(mL3)*(4*pow2(mL3)*(pow2(mL3)*pow2(Mu) - pow4(mL3) +
        pow4(Mu))*pow6(mE3) + 18*pow6(mL3)*pow6(Mu) + 2*(-2*pow2(mL3)*pow2(Mu)
        + pow4(mL3))*pow8(mE3) - 3*pow4(Mu)*pow8(mL3) - 25*pow4(mL3)*pow8(Mu) -
        pow4(mE3)*(41*pow4(mL3)*pow4(Mu) - 28*pow2(Mu)*pow6(mL3) - 18*pow2(mL3)
        *pow6(Mu) + 4*pow8(mL3) + pow8(Mu)) + 4*pow2(mE3)*(-11*pow4(Mu)*pow6(
        mL3) + 15*pow4(mL3)*pow6(Mu) + 2*pow2(Mu)*pow8(mL3) - 7*pow2(mL3)*pow8(
        Mu)) + 12*pow2(mL3)*power10(Mu))))/(pow2(mE3)*pow2(mL3)*pow2(pow2(mE3)
        - pow2(Mu))*pow2(pow2(mL3) - pow2(Mu)))))/pow4(pow2(mE3) - pow2(mL3)))/
        pow2(cbeta));
}

double ThresholdCalculator::getDeltaLambdaYt2Yb4(int omitLogs) const
{
   using std::log;
   using himalaya::dilog;

   const double MR2 = pow2(p.scale);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mU32 = p.mu2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mU3 = std::sqrt(mU32);
   const double mD3 = std::sqrt(mD32);
   const double mA = p.MA;
   const double beta = std::atan(p.vu/p.vd);
   const double sbeta = std::sin(beta);
   const double cbeta = std::cos(beta);
   const double Xt = p.Au(2,2) - p.mu*p.vd/p.vu;
   const double Yt = p.Au(2,2) + p.mu*p.vu/p.vd;
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double Yb = p.Ad(2,2) + p.mu*p.vd/p.vu;
   const double Mu = p.mu;
   const double Mu2 = pow2(Mu);
   const double Xt2 = pow2(Xt);
   const double mA2 = pow2(mA);
   const double Xb2 = pow2(Xb);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);
   const double lmUMR = omitLogs*log(Mu2 / MR2);
   const double eps = mQ3*0.01;

   if(std::abs(mU3 - mQ3) < eps && std::abs(mD3 - mQ3) < eps && std::abs(Mu - mQ3) < eps){
     const double M = std::sqrt(mQ3*mD3);
     return (2.675765602178715*pow2(cbeta)*(cbeta*(1.4948992530373513*sbeta*pow3(M)*
        pow3(Xb) - 8.96939551822411*sbeta*Xt*pow5(M)) + pow4(cbeta)*(Xb2*pow2(
        M)*(-0.3737248132593378*Xb2 - 0.34273774426277803*Xb*Yt - 
        0.514106616394167*pow2(Yt)) + (3.5039051224688698*Xb2 + 
        5.046224971651371*Xb*Yt + 1.2615562429128429*pow2(Yt))*pow4(M) + 
        0.05195844587770295*pow2(Yt)*pow4(Xb) + 1.*pow6(M)) + pow2(sbeta)*(-
        0.20235594112794902*Xb2*Xt2*pow2(sbeta)*pow2(Yb) + pow2(M)*(
        1.4948992530373513*Xb*Xt2*Yb*pow2(sbeta) + 0.7474496265186756*Xt2*pow2(
        sbeta)*pow2(Yb) + Xb2*(1.1211744397780137*Xt2 + 1.4948992530373513*Xt*
        Yb*pow2(sbeta) + 0.7474496265186756*pow2(sbeta)*pow2(Yb))) + (Xb*(-
        7.0078102449377395*Xt - 7.0078102449377395*Yb)*pow2(sbeta) - 
        7.0078102449377395*Xt*Yb*pow2(sbeta) + Xb2*(-2.2423488795560274 + 
        1.2615562429128426*pow2(sbeta)) + Xt2*(-2.2423488795560274 + 
        1.2615562429128426*pow2(sbeta)) - 3.5039051224688698*pow2(sbeta)*pow2(
        Yb))*pow4(M) + (-17.93879103644822 - 10.511715367406609*pow2(sbeta))*
        pow6(M)) + pow2(cbeta)*(Xt*(0.7474496265186756*Xb*Xt - 
        0.6125456657667098*Yb*Yt)*pow2(sbeta)*pow3(Xb) + Xb*pow2(M)*(Xb2*(-
        0.342737744262778*Yb - 0.34273774426277803*Yt)*pow2(sbeta) + 
        2.9897985060747025*Xt*Yb*Yt*pow2(sbeta) + Xb*Xt*(-6.727046638668082*Xt
        + 1.4948992530373513*Yb + 1.4948992530373513*Yt)*pow2(sbeta) + (
        0.3737248132593378 + 0.3737248132593378*pow2(sbeta))*pow3(Xb)) + (-
        8.407868305684794*Xb2 + Xt*(8.96939551822411*Xt - 7.0078102449377395*Yb
        - 7.0078102449377395*Yt) + Xb*(-7.0078102449377395*Xt + 
        5.0462249716513705*Yb + 5.0462249716513705*Yt))*pow2(sbeta)*pow4(M) + (
        8.969395518224111 + 6.242348879556027*pow2(sbeta))*pow6(M))))/pow6(M);
   }

   return -6 + (48*Xb*std::abs(Mu)*(dilog(1 - Mu2/mD32) - dilog(1 - Mu2/mQ32)))/(cbeta*(
        mD32 - mQ32)*sbeta) + (48*Xt*std::abs(Mu)*((mQ32 + Mu2)*(Mu2 - mU32)*dilog(1
        - Mu2/mQ32) + (mQ32 - Mu2)*(Mu2 + mU32)*dilog(1 - Mu2/mU32)))/(cbeta*(-
        mQ32 + Mu2)*(mQ32 - mU32)*(Mu2 - mU32)*sbeta) - (24*dilog(1 - Mu2/mD32)
        )/pow2(cbeta) + 8*pow2(Pi) + 12*dilog(1 - Mu2/mQ32)*(-2/pow2(cbeta) +
        1/pow2(sbeta)) - 3/pow2(sbeta) + (12*Mu2)/(mQ32*pow2(sbeta)) - (6*mU32)
        /(mQ32*pow2(sbeta)) + (6*Mu2)/((-Mu2 + mU32)*pow2(sbeta)) + (12*dilog(1
        - Mu2/mU32))/pow2(sbeta) + (3*pow2(cbeta))/pow2(sbeta) - (6*mA2*pow2(
        cbeta))/(mQ32*pow2(sbeta)) + (2*pow2(cbeta)*pow2(Pi))/pow2(sbeta) - (6*
        pow2(cbeta)*pow2(Yt))/(mQ32*pow2(sbeta)) + (24 + (6*pow2(cbeta))/pow2(
        sbeta))*pow2(log(mA2/MR2)) + (48*((-mD32 - mQ32 + 2*Mu2)*sbeta*std::abs(Mu)*
        dilog(1 - Mu2/mD32) + (mD32 + mQ32 - 2*Mu2)*sbeta*std::abs(Mu)*dilog(1 -
        Mu2/mQ32) - 2*(mD32 - mQ32)*(2*sbeta*std::abs(Mu) + cbeta*(Yb + Yt)*pow2(
        sbeta) + Yt*pow3(cbeta)))*pow3(Xb))/(cbeta*pow2(sbeta)*pow3(mD32 -
        mQ32)) - (24*dilog(1 - Mu2/mU32)*pow4(Mu))/(pow2(Mu2 - mU32)*pow2(
        sbeta)) + pow2(log(mU32/mQ32))*((24*(Mu2 + mU32)*Xt*std::abs(Mu))/(cbeta*(
        Mu2 - mU32)*(-mQ32 + mU32)*sbeta) - (6*Xt2*(2*mU32 + 3*(mQ32 - mU32 +
        Xb2)*pow2(cbeta)))/(pow2(cbeta)*pow2(mQ32 - mU32)) + (6 - (12*pow4(Mu))
        /pow2(Mu2 - mU32))/pow2(sbeta)) + 6*pow2(lmQ3MR)*(3 + (8*pow(Mu2,1.5)
        *Xt)/(cbeta*(-mQ32 + Mu2)*(Mu2 - mU32)*sbeta) + (2*(-1 + pow2(sbeta)))/
        pow2(cbeta) + pow2(cbeta)/pow2(sbeta) + (2 - (2*pow4(Mu))/pow2(Mu2 -
        mU32))/pow2(sbeta)) + (12*Xb2*(3*mD32*mQ32*pow2(cbeta) - 2*mD32*Mu2*
        pow2(cbeta) + 2*mQ32*Mu2*pow2(cbeta) + mD32*mU32*pow2(cbeta) - mQ32*
        mU32*pow2(cbeta) - 2*mD32*mQ32*dilog(1 - Mu2/mQ32)*pow2(cbeta) + 2*
        mQ32*Mu2*dilog(1 - Mu2/mQ32)*pow2(cbeta) + 2*mD32*mQ32*pow2(sbeta) - 4*
        mD32*mQ32*pow2(cbeta)*pow2(sbeta) + 6*(mD32 - mQ32)*mQ32*dilog(1 -
        mD32/mQ32)*pow2(cbeta)*pow2(sbeta) + 2*mQ32*dilog(1 - Mu2/mD32)*((mD32
        - Mu2)*pow2(cbeta) + mQ32*pow2(sbeta)) + mA2*mD32*pow4(cbeta) - mA2*
        mQ32*pow4(cbeta) - 2*mD32*mQ32*pow4(cbeta) + mD32*pow2(Yt)*pow4(cbeta)
        - mQ32*pow2(Yt)*pow4(cbeta) - 3*pow2(cbeta)*pow4(mQ3) - 2*pow2(sbeta)*
        pow4(mQ3) - 2*dilog(1 - Mu2/mQ32)*pow2(sbeta)*pow4(mQ3) + 4*pow2(cbeta)
        *pow2(sbeta)*pow4(mQ3) + 2*pow4(cbeta)*pow4(mQ3) - 2*mD32*mQ32*pow4(
        sbeta) + 2*pow4(mQ3)*pow4(sbeta)))/(mQ32*pow2(cbeta)*pow2(mD32 - mQ32)*
        pow2(sbeta)) + 6*lmUMR*Mu2*((4*(1/(mD32 - Mu2) + 1/(mQ32 - Mu2)))/pow2(
        cbeta) + (4*Xt2)/((mQ32 - Mu2)*(mQ32 - mU32)*pow2(cbeta)) + (4*Xb2*(1/(
        (mD32 - Mu2)*pow2(cbeta)) + 1/(mQ32*pow2(sbeta))))/(mD32 - mQ32) - (2/
        mQ32 + (Mu2 + 2*mU32)/pow2(Mu2 - mU32))/pow2(sbeta) + (2*(mD32 + 2*
        mQ32)*pow4(Xb))/(mQ32*pow2(sbeta)*pow3(-mD32 + mQ32))) + 6*pow2(log(
        mD32/mQ32))*((4*Xb*std::abs(Mu))/(cbeta*mD32*sbeta - cbeta*mQ32*sbeta) - 2/
        pow2(cbeta) + (2*Xb2*((mD32 - Mu2)*pow2(cbeta) + mQ32*pow2(sbeta)))/(
        pow2(cbeta)*pow2(mD32 - mQ32)*pow2(sbeta)) - (4*(mD32 + mQ32 - 2*Mu2)*
        std::abs(Mu)*pow3(Xb))/(cbeta*sbeta*pow3(mD32 - mQ32)) - ((mD32 + 2*mQ32 -
        3*Mu2)*(mD32 - Mu2)*pow4(Xb))/(pow2(sbeta)*pow4(mD32 - mQ32))) + (6*
        phixyz(mA2,mU32,mD32)*((2*Xb*Yt*(pow2(cbeta) + pow2(sbeta))*(-3*mD32*
        mU32 - mA2*(3*mD32 + 2*mU32) + pow4(mA) + 2*pow4(mD3) + pow4(mU3)))/((
        mD32 - mQ32)*pow2(sbeta)) - (2*(mD32 - mQ32 + Xb2)*Xt*Yb*((mD32 - mQ32
        + Xb*Yt)*pow2(cbeta) + (mD32 - mQ32)*pow2(sbeta))*(-3*mD32*mU32 - mA2*(
        3*mD32 + 2*mU32) + pow4(mA) + 2*pow4(mD3) + pow4(mU3)))/((mQ32 - mU32)*
        pow2(cbeta)*pow2(mD32 - mQ32)) + ((mD32 - mQ32 + Xb2)*Xt2*pow2(sbeta)*
        pow2(Yb)*(-3*mD32*mU32 - mA2*(3*mD32 + 2*mU32) + pow4(mA) + 2*pow4(mD3)
        + pow4(mU3)))/((mD32 - mQ32)*pow2(cbeta)*pow2(mQ32 - mU32)) + (pow2(
        pow2(cbeta) + pow2(sbeta))*(-3*mD32*mU32 - mA2*(3*mD32 + 2*mU32) +
        pow4(mA) + 2*pow4(mD3) + pow4(mU3)))/(pow2(cbeta)*pow2(sbeta)) + (2*Yt*
        (pow2(cbeta) + pow2(sbeta))*pow3(Xb)*(-3*mD32*mU32 - mA2*(3*mD32 + 2*
        mU32) + pow4(mA) + 2*pow4(mD3) + pow4(mU3)))/(pow2(mD32 - mQ32)*pow2(
        sbeta)) + (Xb2*(-3*mD32*mU32 - mA2*(3*mD32 + 2*mU32) + pow4(mA) + 2*
        pow4(mD3) + pow4(mU3))*(2*(mD32 - mQ32)*pow2(cbeta)*pow2(sbeta) + (mD32
        - mQ32 + pow2(Yt))*pow4(cbeta) + (mD32 - mQ32)*pow4(sbeta)))/(pow2(
        cbeta)*pow2(mD32 - mQ32)*pow2(sbeta)) + (pow2(cbeta)*pow2(Yt)*(-3*mD32*
        mU32 - mA2*(3*mD32 + 2*mU32) + pow4(mA) + 2*pow4(mD3) + pow4(mU3))*
        pow4(Xb))/(pow2(sbeta)*pow3(mD32 - mQ32)) + deltaxyz(mA2,mU32,mD32)*((-
        2*Xb*Yt*(pow2(cbeta) + pow2(sbeta)))/((mD32 - mQ32)*pow2(sbeta)) -
        pow2(pow2(cbeta) + pow2(sbeta))/(pow2(cbeta)*pow2(sbeta)) - (2*(3*mD32
        - mQ32)*Yt*(pow2(cbeta) + pow2(sbeta))*pow3(Xb))/(pow2(sbeta)*pow3(mD32
        - mQ32)) + (2*Xt*(-(Xb*pow2(mD32 - mQ32)*((mD32 - Yb*Yt)*pow2(cbeta) +
        mD32*pow2(sbeta))) + (mD32 - mQ32)*Xb2*(-((mQ32*Yb + mD32*(-2*Yb + Yt))
        *pow2(cbeta)) + (2*mD32 - mQ32)*Yb*pow2(sbeta)) + Yb*(pow2(cbeta) +
        pow2(sbeta))*pow3(mD32 - mQ32) + (3*mD32 - mQ32)*Yb*Yt*pow2(cbeta)*
        pow3(Xb)))/((mQ32 - mU32)*pow2(cbeta)*pow3(mD32 - mQ32)) - (Xt2*Yb*
        pow2(sbeta)*(mQ32*(mQ32 - Xb2)*Yb + 2*mD32*(mQ32*(Xb - Yb) + Xb2*Yb) +
        (-2*Xb + Yb)*pow4(mD3)))/(pow2(cbeta)*pow2(mD32 - mQ32)*pow2(mQ32 -
        mU32)) - (Xb2*(2*(2*mD32 - mQ32)*pow2(cbeta)*pow2(sbeta) + (2*mD32 -
        mQ32 + pow2(Yt))*pow4(cbeta) + (2*mD32 - mQ32)*pow4(sbeta)))/(pow2(
        cbeta)*pow2(mD32 - mQ32)*pow2(sbeta)) - ((4*mD32 - mQ32)*pow2(cbeta)*
        pow2(Yt)*pow4(Xb))/(pow2(sbeta)*pow4(mD32 - mQ32)))))/pow4(mD3) + (6*
        phixyz(mA2,mQ32,mD32)*(((mA2 - mQ32)*pow2(sbeta)*pow2(Yb))/pow2(cbeta)
        + (Xb2*(mA2*(2*mD32 - mQ32) + pow2(mD32 - mQ32))*pow2(sbeta)*pow2(Yb))/
        (pow2(cbeta)*pow2(mD32 - mQ32)) + (deltaxyz(mA2,mQ32,mD32)*(2*mD32*(
        mD32 - mQ32)*Xb2 - 2*Xb*Yb*pow2(mD32 - mQ32) - 2*(3*mD32 - mQ32)*Yb*
        pow3(Xb) - (2*Xt*(-(Xb*pow2(mD32 - mQ32)*((mD32 - Yb*Yt)*pow2(cbeta) +
        mD32*pow2(sbeta))) - (mD32 - mQ32)*Xb2*((mQ32*Yb + mD32*(-2*Yb + Yt))*
        pow2(cbeta) + (-2*mD32 + mQ32)*Yb*pow2(sbeta)) + Yb*(pow2(cbeta) +
        pow2(sbeta))*pow3(mD32 - mQ32) + (3*mD32 - mQ32)*Yb*Yt*pow2(cbeta)*
        pow3(Xb)))/((mQ32 - mU32)*pow2(cbeta)) + ((mD32 - mQ32)*Xt2*Yb*pow2(
        sbeta)*(mQ32*(mQ32 - Xb2)*Yb + 2*mD32*(mQ32*(Xb - Yb) + Xb2*Yb) + (-2*
        Xb + Yb)*pow4(mD3)))/(pow2(cbeta)*pow2(mQ32 - mU32))))/pow3(mD32 -
        mQ32) + (2*(mD32 - mQ32 + Xb2)*Xt*Yb*((mD32 - mQ32 + Xb*Yt)*pow2(cbeta)
        + (mD32 - mQ32)*pow2(sbeta))*(-3*mD32*mQ32 - mA2*(3*mD32 + 2*mQ32) +
        pow4(mA) + 2*pow4(mD3) + pow4(mQ3)))/((mQ32 - mU32)*pow2(cbeta)*pow2(
        mD32 - mQ32)) + (2*Yb*pow3(Xb)*(-3*mD32*mQ32 - mA2*(3*mD32 + 2*mQ32) +
        pow4(mA) + 2*pow4(mD3) + pow4(mQ3)))/pow2(mD32 - mQ32) + (2*Xb*Yb*(-(
        mD32*(mA2 + mD32 - mQ32)*pow2(sbeta)) + pow2(cbeta)*(-3*mD32*mQ32 -
        mA2*(3*mD32 + 2*mQ32) + pow4(mA) + 2*pow4(mD3) + pow4(mQ3))))/((mD32 -
        mQ32)*pow2(cbeta)) - (Xt2*Yb*pow2(sbeta)*(2*mD32*(mD32 - mQ32)*(mA2 +
        mD32 - mQ32)*(mQ32 - mU32)*Xb + Xb2*Yb*((2*mD32 - 2*mQ32 + mU32)*pow2(
        mD32 - mQ32) + (mD32 - mQ32)*pow4(mA) - mA2*(mD32*(mQ32 - 2*mU32) +
        mQ32*mU32 + 3*pow4(mD3) - 3*pow4(mQ3))) + Yb*pow2(mD32 - mQ32)*(-3*
        mD32*mQ32 - mQ32*mU32 + mA2*(-3*mD32 - 3*mQ32 + mU32) + pow4(mA) + 2*
        pow4(mD3) + 2*pow4(mQ3))))/(pow2(cbeta)*pow2(mD32 - mQ32)*pow2(mQ32 -
        mU32)) + ((mD32 - mQ32 + Xb2)*(mQ32 - mU32 + Xt2)*pow2(sbeta)*pow2(Yb)*
        (pow2(-(mD32*mQ3) + pow3(mQ3)) + 3*mQ32*pow4(mA) + mA2*(2*mD32*mQ32 +
        pow4(mD3) - 3*pow4(mQ3)) - pow6(mA)))/((mD32 - mQ32)*(mQ32 - mU32)*
        deltaxyz(mA2,mQ32,mD32)*pow2(cbeta))))/pow4(mD3) + (3*log(mA2/MR2)*((
        pow4(cbeta)*(2*mQ32*pow2(Yt) - mU32*pow2(Yt) + mA2*(2*mQ32 + pow2(Yt))
        + 2*pow4(mQ3)))/pow4(mQ3) - (2*(mD32 - mQ32 + Xb2)*Xt2*pow2(Yb)*pow4(
        sbeta))/(mD32*(mD32 - mQ32)*(mQ32 - mU32)) + (2*(mD32 - mQ32 + Xb2)*(
        mQ32 - mU32 + Xt2)*pow2(mA2 + mD32 - mQ32)*pow2(Yb)*pow4(sbeta))/(mD32*
        (mD32 - mQ32)*(mQ32 - mU32)*deltaxyz(mA2,mQ32,mD32)) + 2*(4 - pow2(Yb)/
        mD32)*pow4(sbeta) - (2*Xb2*(mD32*((2*mQ32 - mU32)*pow2(Yt) + mA2*(2*
        mQ32 + pow2(Yt)))*pow4(cbeta) + pow2(Yb)*pow4(mQ3)*pow4(sbeta)))/(mD32*
        (mD32 - mQ32)*pow4(mQ3)) + (pow4(cbeta)*((mD32 - mQ32)*(2*mQ32 - mU32)*
        pow2(Yt) + mA2*(-(mQ32*pow2(Yt)) + mD32*(2*mQ32 + pow2(Yt)) + 10*pow4(
        mQ3)))*pow4(Xb))/(pow3(mD32 - mQ32)*pow4(mQ3)) - (pow2(-mD32 + mQ32 +
        Xb2)*pow2(Yt)*pow4(cbeta)*((2*mQ32 - mU32)*pow2(mQ32 - mU32) - (2*mQ32
        + 3*mU32)*pow4(mA) - mA2*(2*mQ32*mU32 + pow4(mQ3) - 3*pow4(mU3)) +
        pow6(mA)))/(deltaxyz(mA2,mU32,mQ32)*pow2(mD32 - mQ32)*pow4(mQ3))))/(
        pow2(cbeta)*pow2(sbeta)) - (6*pow4(Xb)*(-14*mD32*mQ32*Mu2 + 4*mD32*
        mQ32*mU32 + 2*mQ32*(mD32 + 2*mQ32 - 3*Mu2)*(mD32 - Mu2)*dilog(1 - Mu2/
        mD32) + 8*mD32*mQ32*Mu2*dilog(1 - Mu2/mQ32) + 4*mA2*mD32*mQ32*pow2(
        cbeta) + 6*mQ32*dilog(1 - mD32/mQ32)*pow2(mD32 - mQ32)*pow2(sbeta) +
        10*mD32*mQ32*pow2(cbeta)*pow2(Yt) + 9*mQ32*pow4(mD3) - 2*Mu2*pow4(mD3)
        + mU32*pow4(mD3) - 2*mQ32*dilog(1 - Mu2/mQ32)*pow4(mD3) + mA2*pow2(
        cbeta)*pow4(mD3) - 16*mQ32*pow2(sbeta)*pow4(mD3) + pow2(cbeta)*pow2(Yt)
        *pow4(mD3) - 6*mD32*pow4(mQ3) + 16*Mu2*pow4(mQ3) - 5*mU32*pow4(mQ3) -
        4*mD32*dilog(1 - Mu2/mQ32)*pow4(mQ3) + 4*Mu2*dilog(1 - Mu2/mQ32)*pow4(
        mQ3) - 5*mA2*pow2(cbeta)*pow4(mQ3) + 32*mD32*pow2(sbeta)*pow4(mQ3) -
        11*pow2(cbeta)*pow2(Yt)*pow4(mQ3) - 6*mQ32*dilog(1 - Mu2/mQ32)*pow4(Mu)
        - 3*pow6(mQ3) - 16*pow2(sbeta)*pow6(mQ3)))/(mQ32*pow2(sbeta)*pow4(mD32
        - mQ32)) + (6*phixyz(mA2,mQ32,mQ32)*((2*mA2*(mA2 - 5*mQ32)*(-mD32 +
        mQ32 + Xb2)*Xt*(-mD32 + mQ32 + Xb*Yb)*Yt)/((mQ32 - mU32)*pow2(mD32 -
        mQ32)) - (mA2*mQ32*(10*pow2(cbeta) + pow2(sbeta)))/pow2(cbeta) + (2*
        mA2*Xb*Yb*((mA2 - 5*mQ32)*pow2(cbeta) - mQ32*pow2(sbeta)))/((-mD32 +
        mQ32)*pow2(cbeta)) - (mA2*Xb2*(2*(mA2 - 5*mQ32)*(mD32 - mQ32)*pow2(
        cbeta) + mQ32*pow2(sbeta)*pow2(Yb)))/(pow2(cbeta)*pow2(mD32 - mQ32)) -
        (mA2*mQ32*Xt2*pow2(sbeta)*pow2(-mD32 + mQ32 + Xb*Yb))/((mQ32 - mU32)*
        pow2(cbeta)*pow2(mD32 - mQ32)) + (2*mA2*(mA2 - 5*mQ32)*Yb*pow3(Xb))/
        pow2(mD32 - mQ32) + 2*pow4(mA) + deltaxyz(mA2,mQ32,mQ32)*(-2 + (2*Xb*
        Yb)/(mD32 - mQ32) + (2*(mD32 - 2*mQ32)*Xb2)/pow2(mD32 - mQ32) - (mQ32*
        Xt2*pow2(sbeta)*pow2(-mD32 + mQ32 + Xb*Yb))/(pow2(cbeta)*pow2(mD32 -
        mQ32)*pow2(mQ32 - mU32)) + (2*(mD32 - 3*mQ32)*Yb*pow3(Xb))/pow3(-mD32 +
        mQ32) + (2*Xt*(mQ32*(-mD32 + mQ32)*Xb*(-mD32 + mQ32 + Xb*Yb)*pow2(
        sbeta) + pow2(cbeta)*(-3*mQ32*Yb*Yt*pow3(Xb) + (mQ32*(Xb - 3*Yt) - Xb*(
        Xb + Yb)*Yt)*pow4(mD3) + Xb*(Xb*Yb - 2*Xb*Yt - Yb*Yt)*pow4(mQ3) + mD32*
        (mQ32*Xb*(-(Xb*Yb) + 3*Xb*Yt + 2*Yb*Yt) + Yb*Yt*pow3(Xb) + (-2*Xb + 3*
        Yt)*pow4(mQ3)) + Yt*pow6(mD3) + (Xb - Yt)*pow6(mQ3))))/((mQ32 - mU32)*
        pow2(cbeta)*pow3(-mD32 + mQ32)))))/pow4(mQ3) + lmQ3MR*((-6*(mD32 - mQ32
        + Xb2)*(mQ32 - mU32 + Xt2)*pow2(mA2 + mD32 - mQ32)*pow2(sbeta)*pow2(Yb)
        )/(mD32*(mD32 - mQ32)*(mQ32 - mU32)*deltaxyz(mA2,mQ32,mD32)*pow2(cbeta)
        ) + (48*(2*sbeta*std::abs(Mu) + cbeta*(Yb + Yt)*pow2(sbeta) + Yt*pow3(cbeta)
        )*pow3(Xb))/(cbeta*pow2(mD32 - mQ32)*pow2(sbeta)) + (6*Xb2*(4 + (-4 - (
        2*mU32)/mQ32)/pow2(sbeta) + (-(Mu2*pow2(sbeta)*pow2(Yb)) + mD32*(-2*
        Mu2*(1 + pow2(sbeta)) + pow2(sbeta)*pow2(Yb)) + 2*(-1 + pow2(sbeta))*
        pow4(mD3))/(mD32*(mD32 - Mu2)*pow2(cbeta)) + (pow2(cbeta)*((mA2 - mU32)
        *pow2(Yt) + 2*pow4(mQ3)))/(pow2(sbeta)*pow4(mQ3))))/(mD32 - mQ32) + (
        12*lmUMR*((-4*pow(Mu2,1.5)*(Mu2 - mU32)*sbeta*Xt)/(cbeta*(-mQ32 +
        Mu2)) + pow4(Mu)))/(pow2(Mu2 - mU32)*pow2(sbeta)) + 3*((8*Mu2*(1/(-mD32
        + Mu2) + 1/(-mQ32 + Mu2)) + 2*pow2(sbeta)*(-4 + pow2(Yb)/mD32))/pow2(
        cbeta) - (pow2(cbeta)*((mA2 - mU32)*pow2(Yt) + 4*pow4(mQ3)))/(pow2(
        sbeta)*pow4(mQ3)) + (2*pow2(-(Mu2*mU3) + pow3(mU3)) - 2*mQ32*(-4*Mu2*
        mU32*(1 + 3*pow2(sbeta)) + (1 + 6*pow2(sbeta))*pow4(Mu) + 6*pow2(sbeta)
        *pow4(mU3)))/(mQ32*pow2(Mu2 - mU32)*pow2(sbeta))) + (3*pow2(cbeta)*
        pow2(-mD32 + mQ32 + Xb2)*pow2(Yt)*((2*mQ32 - mU32)*pow2(mQ32 - mU32) -
        (2*mQ32 + 3*mU32)*pow4(mA) - mA2*(2*mQ32*mU32 + pow4(mQ3) - 3*pow4(mU3)
        ) + pow6(mA)))/(deltaxyz(mA2,mU32,mQ32)*pow2(mD32 - mQ32)*pow2(sbeta)*
        pow4(mQ3)) - (3*pow4(Xb)*(mQ32*(mA2 - mU32)*pow2(cbeta)*pow2(Yt) - 2*(
        6*Mu2 - 5*mU32 - 6*pow2(cbeta)*pow2(Yt))*pow4(mQ3) + mD32*(2*mQ32*mU32
        + (-mA2 + mU32)*pow2(cbeta)*pow2(Yt) - 2*(-5 + 8*pow2(sbeta))*pow4(mQ3)
        ) + 2*(1 + 8*pow2(sbeta))*pow6(mQ3)))/(pow2(sbeta)*pow3(-mD32 + mQ32)*
        pow4(mQ3)) + Xt2*((6*Xb2*(2*mD32*(-mQ32 + mU32)*pow2(cbeta) + mQ32*
        pow2(sbeta)*pow2(Yb)))/(mD32*(mD32 - mQ32)*mQ32*(mQ32 - mU32)*pow2(
        cbeta)) - (6*(mD32 + 5*mQ32)*pow4(Xb))/(mQ32*pow3(-mD32 + mQ32)) + (6*(
        ((mQ32 - mU32)*pow2(sbeta)*pow2(Yb))/mD32 - ((2*mU32*(-1 + 5*pow2(
        cbeta) + pow2(sbeta)) + Mu2*(2 - 3*pow2(cbeta) + 2*pow2(sbeta)))*pow4(
        mQ3) + Mu2*pow2(cbeta)*pow4(mU3) - mQ32*(2*Mu2*mU32*(1 + 5*pow2(cbeta)
        + pow2(sbeta)) + pow2(cbeta)*pow4(mU3)) + (2 + 3*pow2(cbeta) - 2*pow2(
        sbeta))*pow6(mQ3))/(mQ32*(mQ32 - Mu2))))/(pow2(cbeta)*pow2(mQ32 - mU32)
        ))) + log(mU32/mQ32)*((6*Xb2*((-mA2 + mU32)*pow2(cbeta)*pow2(Yt) + 2*
        mQ32*(mU32 + pow2(cbeta)*pow2(Yt))))/((-mD32 + mQ32)*pow2(sbeta)*pow4(
        mQ3)) + (12*lmUMR*((-4*pow(Mu2,1.5)*(Mu2 - mU32)*sbeta*Xt)/(cbeta*(-
        mQ32 + mU32)) + pow4(Mu)))/(pow2(Mu2 - mU32)*pow2(sbeta)) + (3*((2*(
        mU32 + pow2(cbeta)*pow2(Yt)))/mQ32 + ((-mA2 + mU32)*pow2(cbeta)*pow2(
        Yt))/pow4(mQ3) + (2*(2*Mu2*mU32 + pow4(mU3)))/pow2(Mu2 - mU32)))/pow2(
        sbeta) + 6*lmQ3MR*(2 + pow2(sbeta)/pow2(cbeta) + (2*Xt*(4*cbeta*(Mu2 +
        mU32)*std::abs(Mu) + (Mu2 - mU32)*sbeta*(Xb + Yb + Yt)*pow2(cbeta) + (Mu2 -
        mU32)*(Xb + Yb)*pow3(sbeta)))/((mQ32 - mU32)*(-Mu2 + mU32)*sbeta*pow2(
        cbeta)) + (2 + pow2(cbeta) - (4*pow4(Mu))/pow2(Mu2 - mU32))/pow2(sbeta)
        - (Xt2*(mQ32*(2*mU32*(1 + 5*pow2(cbeta)) + (mA2 - Yb*(2*Xb + Yb))*pow2(
        sbeta)) + mU32*pow2(sbeta)*(-mA2 + 2*Xb*Yb + pow2(Yb)) - pow2(sbeta)*
        pow4(mQ3) + (-2 + 2*pow2(cbeta) + pow2(sbeta))*pow4(mU3)))/(pow2(cbeta)
        *pow3(mQ32 - mU32))) - (6*log(mA2/MR2)*(2*(mQ32 - mU32)*(mQ32 - mU32 -
        Xt*(Xb + Yb + Yt))*pow2(cbeta)*pow2(sbeta) + pow2(mQ32 - mU32)*pow4(
        cbeta) + (mQ32*(-2*mU32 + Xt*(-2*Xb + Xt - 2*Yb)) + mU32*Xt*(2*Xb - Xt
        + 2*Yb) + Xt2*(-mA2 + 2*Xb*Yb + pow2(Yb)) + pow4(mQ3) + pow4(mU3))*
        pow4(sbeta)))/(pow2(cbeta)*pow2(mQ32 - mU32)*pow2(sbeta)) - (3*(mQ32*(
        mA2 - mU32)*pow2(cbeta)*pow2(Yt) + mD32*((-mA2 + mU32)*pow2(cbeta)*
        pow2(Yt) + 2*mQ32*(mU32 + pow2(cbeta)*pow2(Yt))) + 2*(5*mU32 - pow2(
        cbeta)*pow2(Yt))*pow4(mQ3))*pow4(Xb))/(pow2(sbeta)*pow3(-mD32 + mQ32)*
        pow4(mQ3)) + (Xt2*(((mQ32 - mU32)*(-6*mQ32*mU32*(-4 + 5*pow2(cbeta) +
        4*pow2(sbeta)) + 6*pow2(cbeta)*pow4(mU3)))/(mQ32*pow2(cbeta)) + ((mQ32
        - mU32)*Xb2*(-12*mQ32*mU32*(1 + pow2(cbeta)) + 12*pow2(cbeta)*pow4(mU3)
        ))/(mQ32*(-mD32 + mQ32)*pow2(cbeta)) - (6*mU32*(-2*mQ32*mU32 + 13*pow4(
        mQ3) + pow4(mU3))*pow4(Xb))/pow2(-(mD32*mQ3) + pow3(mQ3))))/pow3(mQ32 -
        mU32) + (3*pow2(cbeta)*pow2(-mD32 + mQ32 + Xb2)*pow2(Yt)*(-((4*mQ32 +
        3*mU32)*pow4(mA)) + mU32*pow4(mQ3) + 2*mQ32*pow4(mU3) + mA2*(2*mQ32*
        mU32 + 5*pow4(mQ3) + 3*pow4(mU3)) + pow6(mA) - 2*pow6(mQ3) - pow6(mU3))
        )/(deltaxyz(mA2,mU32,mQ32)*pow2(mD32 - mQ32)*pow2(sbeta)*pow4(mQ3))) +
        log(mD32/mQ32)*(12*(-2 + ((2*Mu2)/(-mD32 + Mu2) - pow2(sbeta))/pow2(
        cbeta) - pow2(cbeta)/pow2(sbeta)) - (12*(mA2 + mD32 - mQ32)*(mD32 -
        mQ32 + Xb2)*(mQ32 - mU32 + Xt2)*pow2(sbeta)*pow2(Yb))/((mD32 - mQ32)*(
        mQ32 - mU32)*deltaxyz(mA2,mQ32,mD32)*pow2(cbeta)) - (24*Xb*(2*sbeta*
        std::abs(Mu) + cbeta*(Yb + Yt)*pow2(sbeta) + Yt*pow3(cbeta)))/(cbeta*(mD32 -
        mQ32)*pow2(sbeta)) + (24*(3*mD32 + mQ32)*(2*sbeta*std::abs(Mu) + cbeta*(Yb +
        Yt)*pow2(sbeta) + Yt*pow3(cbeta))*pow3(Xb))/(cbeta*pow2(sbeta)*pow3(
        mD32 - mQ32)) - (12*Xb2*(-((mD32 - Mu2)*pow2(cbeta)*(2*Mu2 - mU32 + 3*
        mD32*(-1 + pow2(sbeta)) + mQ32*pow2(sbeta))) + (mD32 - Mu2)*(mA2 - mD32
        - mQ32 + pow2(Yt))*pow4(cbeta) + pow2(sbeta)*(mQ32*Mu2*(-2 + pow2(
        sbeta)) + mD32*(-mQ32 + Mu2)*pow2(sbeta) - (-2 + pow2(sbeta))*pow4(mD3)
        )))/((mD32 - Mu2)*pow2(cbeta)*pow2(mD32 - mQ32)*pow2(sbeta)) + (12*
        lmUMR*pow3(Xb)*(-8*mD32*pow(Mu2,1.5)*sbeta + 8*mQ32*pow(Mu2,1.5)*
        sbeta + 3*cbeta*Xb*pow4(Mu)))/(cbeta*pow2(sbeta)*pow4(mD32 - mQ32)) + (
        12*(mQ32*(-2*Mu2 + mU32 + mQ32*pow2(sbeta)) + mD32*(-7*Mu2 + 2*mU32 +
        mQ32*(3 + 6*pow2(sbeta))) + pow2(cbeta)*(mA2*(2*mD32 + mQ32) + (5*mD32
        + mQ32)*pow2(Yt)) + (3 - 7*pow2(sbeta))*pow4(mD3))*pow4(Xb))/(pow2(
        sbeta)*pow4(mD32 - mQ32)) - (6*log(mA2/MR2)*(2*(mD32 - mQ32)*Xb*Xt*
        pow2(sbeta)*(-2*Xb2*Yb*Yt*pow2(cbeta) + pow2(mD32 - mQ32)*(pow2(cbeta)
        + pow2(sbeta)) - (mD32 - mQ32)*Xb*((Yb - Yt)*pow2(cbeta) + Yb*pow2(
        sbeta))) - 2*(mD32 - mQ32)*pow2(cbeta)*(-(mQ32*Yt*pow2(cbeta)) + 2*
        mU32*Yt*pow2(cbeta) + mQ32*Yb*pow2(sbeta) - mQ32*Yt*pow2(sbeta) + 2*
        mU32*Yt*pow2(sbeta) - 2*mA2*(Yt*pow2(cbeta) + (Yb + Yt)*pow2(sbeta)) -
        mD32*(Yt*pow2(cbeta) + (Yb + Yt)*pow2(sbeta)))*pow3(Xb) - pow2(pow2(
        cbeta) + pow2(sbeta))*pow4(mD32 - mQ32) + 2*Xb*pow3(mD32 - mQ32)*(-((Yb
        + Yt)*pow2(cbeta)*pow2(sbeta)) - Yt*pow4(cbeta) + Yb*pow4(sbeta)) -
        Xb2*pow2(mD32 - mQ32)*(2*(mD32 - 2*mQ32 + mU32)*pow2(cbeta)*pow2(sbeta)
        + (mA2 - mQ32 + mU32 + pow2(Yt))*pow4(cbeta) + (-mA2 - mQ32 + mU32 +
        pow2(Yb))*pow4(sbeta)) + ((2*mD32 + mQ32 - 3*mU32)*pow2(Yt) + mA2*(4*
        mD32 + 2*mQ32 + 3*pow2(Yt)))*pow4(cbeta)*pow4(Xb)))/(pow2(cbeta)*pow2(
        sbeta)*pow4(mD32 - mQ32)) + (6*lmQ3MR*(2*(mD32 - mQ32)*Xb*Xt*pow2(
        sbeta)*(-2*Xb2*Yb*Yt*pow2(cbeta) + pow2(mD32 - mQ32)*(pow2(cbeta) +
        pow2(sbeta)) - (mD32 - mQ32)*Xb*((Yb - Yt)*pow2(cbeta) + Yb*pow2(sbeta)
        )) + 2*cbeta*(mD32 - mQ32)*(8*pow(Mu2,1.5)*sbeta + 2*cbeta*(-(mU32*
        Yt) + mA2*(Yb + Yt))*pow2(sbeta) + 2*(mA2 - mU32)*Yt*pow3(cbeta) -
        mD32*(4*sbeta*std::abs(Mu) + cbeta*(Yb + Yt)*pow2(sbeta) + Yt*pow3(cbeta)) -
        mQ32*(4*sbeta*std::abs(Mu) + cbeta*(3*Yb + Yt)*pow2(sbeta) + Yt*pow3(cbeta))
        )*pow3(Xb) + pow4(mD32 - mQ32)*(-2*pow2(sbeta) + 2*pow2(cbeta)*pow2(
        sbeta) + pow4(cbeta) + pow4(sbeta)) + 2*Xb*pow3(mD32 - mQ32)*(4*cbeta*
        sbeta*std::abs(Mu) + (Yb + Yt)*pow2(cbeta)*pow2(sbeta) + Yt*pow4(cbeta) +
        Yb*pow4(sbeta)) - Xb2*pow2(mD32 - mQ32)*(mQ32*(-2 + pow2(sbeta))*pow2(
        sbeta) + 2*pow2(cbeta)*(2*Mu2 - mU32 - 2*mQ32*pow2(sbeta) + mU32*pow2(
        sbeta) + mD32*(-2 + 3*pow2(sbeta))) + (-mA2 + mQ32 + mU32 - pow2(Yt))*
        pow4(cbeta) + (-mA2 + mU32 + pow2(Yb))*pow4(sbeta)) - pow2(cbeta)*(-3*
        mA2*pow2(cbeta)*pow2(Yt) + 3*mU32*pow2(cbeta)*pow2(Yt) + mQ32*(-4*Mu2 +
        2*mU32 + pow2(cbeta)*pow2(Yt)) + 2*mD32*(2*mQ32 - 4*Mu2 + 2*mU32 +
        pow2(cbeta)*pow2(Yt)) + (2 - 4*pow2(sbeta))*pow4(mD3) + 4*pow2(sbeta)*
        pow4(mQ3) + 6*pow4(Mu))*pow4(Xb) + 2*Xt2*pow2(cbeta)*pow2(sbeta)*(Xb2*
        pow2(mD32 - mQ32) - (2*mD32 + mQ32)*pow4(Xb))))/(pow2(cbeta)*pow2(
        sbeta)*pow4(mD32 - mQ32)) + (6*Xt2*(-((mQ32*(-2 + 3*pow2(cbeta)) +
        mU32*(2 + 3*pow2(cbeta)))/pow2(cbeta)) - (2*Xb2*(-(mQ32*mU32) + pow4(
        mQ3) + pow2(cbeta)*(-4*mQ32*mU32 + 2*mD32*(2*mQ32 + mU32) - 3*pow4(mQ3)
        + pow4(mU3))))/(pow2(cbeta)*pow2(mD32 - mQ32)) + (pow4(Xb)*((mQ32 + 5*
        mU32)*pow4(mD3) - 5*mU32*pow4(mQ3) + 2*mQ32*pow4(mU3) + 2*mD32*(-12*
        mQ32*mU32 + 7*pow4(mQ3) + 5*pow4(mU3)) - 3*pow6(mQ3)))/pow4(mD32 -
        mQ32)))/pow2(mQ32 - mU32) + log(mU32/mQ32)*((12*Xb*Yt*(pow2(cbeta) +
        pow2(sbeta)))/((mD32 - mQ32)*pow2(sbeta)) + (6*pow2(pow2(cbeta) + pow2(
        sbeta)))/(pow2(cbeta)*pow2(sbeta)) - (12*(-2*mA2 + mD32 + mQ32 + 2*
        mU32)*Yt*(pow2(cbeta) + pow2(sbeta))*pow3(Xb))/(pow2(sbeta)*pow3(mD32 -
        mQ32)) - (12*Xt*(Xb*pow2(mD32 - mQ32)*((-mA2 + mD32 + mU32 + Yb*Yt)*
        pow2(cbeta) + (-mA2 + mD32 + mU32)*pow2(sbeta)) + (mD32 - mQ32)*Xb2*((-
        (mQ32*Yb) - mU32*Yb + mA2*(Yb - Yt) + mD32*Yt + mU32*Yt)*pow2(cbeta) +
        (mA2 - mQ32 - mU32)*Yb*pow2(sbeta)) + Yb*(pow2(cbeta) + pow2(sbeta))*
        pow3(mD32 - mQ32) + (2*mA2 - mD32 - mQ32 - 2*mU32)*Yb*Yt*pow2(cbeta)*
        pow3(Xb)))/((mQ32 - mU32)*pow2(cbeta)*pow3(mD32 - mQ32)) - (6*Xb2*(2*
        pow2(cbeta)*(mU32*(-1 + pow2(sbeta)) + (-mA2 + mQ32)*pow2(sbeta)) + (-
        mA2 + mQ32 + mU32 - pow2(Yt))*pow4(cbeta) + (-mA2 + mQ32 + mU32)*pow4(
        sbeta)))/(pow2(cbeta)*pow2(mD32 - mQ32)*pow2(sbeta)) - (6*(3*(-mA2 +
        mU32)*pow2(cbeta)*pow2(Yt) + 2*mD32*(2*mU32 + pow2(cbeta)*pow2(Yt)) +
        mQ32*(2*mU32 + pow2(cbeta)*pow2(Yt)))*pow4(Xb))/(pow2(sbeta)*pow4(mD32
        - mQ32)) - (Xt2*(-12*(mQ32 - mU32)*(-mA2 + mD32 + mU32)*Xb*Yb*pow2(
        sbeta)*pow3(mD32 - mQ32) + 6*pow4(mD32 - mQ32)*(mU32*pow2(sbeta)*pow2(
        Yb) + mQ32*(mU32*(-2 + 6*pow2(cbeta)) - pow2(sbeta)*pow2(Yb)) + 2*pow4(
        mU3)) - Xb2*pow2(mD32 - mQ32)*(-6*(mQ32 - mU32)*((-mA2 + mU32)*pow2(
        sbeta)*pow2(Yb) + mQ32*(2*mU32 + pow2(sbeta)*pow2(Yb))) + 12*pow2(
        cbeta)*(3*(mQ32 - mU32)*pow4(mD3) + 6*mU32*pow4(mQ3) + mQ32*pow4(mU3) +
        mD32*(-10*mQ32*mU32 + 4*pow4(mU3)) - pow6(mU3))) + 12*pow2(cbeta)*pow4(
        Xb)*(pow4(mQ3)*pow4(mU3) + pow4(mD3)*(-8*mQ32*mU32 + 3*pow4(mQ3) + 2*
        pow4(mU3)) + 3*mU32*pow6(mQ3) - mQ32*pow6(mU3) - 2*mD32*(2*mU32*pow4(
        mQ3) - 3*mQ32*pow4(mU3) + pow6(mU3)))))/(pow2(cbeta)*pow3(mQ32 - mU32)*
        pow4(mD32 - mQ32)))) + (3*phixyz(mA2,mU32,mQ32)*((-4*mQ32*(-mD32 + mQ32
        + Xb2)*Xt*(-mD32 + mQ32 + Xb*Yb)*Yt*(-3*mQ32*mU32 - mA2*(3*mQ32 + 2*
        mU32) + pow4(mA) + 2*pow4(mQ3) + pow4(mU3)))/((mQ32 - mU32)*pow2(mD32 -
        mQ32)) + (4*mQ32*Xb*Yt*(pow2(cbeta) + pow2(sbeta))*(-3*mQ32*mU32 - mA2*
        (3*mQ32 + 2*mU32) + pow4(mA) + 2*pow4(mQ3) + pow4(mU3)))/((-mD32 +
        mQ32)*pow2(sbeta)) + (4*mQ32*Yt*(pow2(cbeta) + pow2(sbeta))*pow3(Xb)*(-
        3*mQ32*mU32 - mA2*(3*mQ32 + 2*mU32) + pow4(mA) + 2*pow4(mQ3) + pow4(
        mU3)))/(pow2(mD32 - mQ32)*pow2(sbeta)) - (pow2(cbeta)*pow2(Yt)*(-6*
        mQ32*mU32 - 6*mA2*(mQ32 + mU32) + 3*pow4(mA) + 2*pow4(mQ3) + 3*pow4(
        mU3)))/pow2(sbeta) + (pow2(cbeta)*pow2(Yt)*pow4(Xb)*((3*mD32 - 7*mQ32)*
        pow4(mA) - 2*mA2*(-7*mQ32*mU32 + 3*mD32*(mQ32 + mU32) - 9*pow4(mQ3)) +
        18*mU32*pow4(mQ3) - 7*mQ32*pow4(mU3) + mD32*(-6*mQ32*mU32 + 2*pow4(mQ3)
        + 3*pow4(mU3)) - 10*pow6(mQ3)))/(pow2(sbeta)*pow3(-mD32 + mQ32)) + (2*
        Xb2*pow2(cbeta)*pow2(Yt)*((3*mD32 - 4*mQ32)*pow4(mA) + 9*mU32*pow4(mQ3)
        + mA2*(8*mQ32*mU32 - 6*mD32*(mQ32 + mU32) + 9*pow4(mQ3)) - 4*mQ32*pow4(
        mU3) + mD32*(-6*mQ32*mU32 + 2*pow4(mQ3) + 3*pow4(mU3)) - 4*pow6(mQ3)))/
        (pow2(mD32 - mQ32)*pow2(sbeta)) + 2*deltaxyz(mA2,mU32,mQ32)*((2*mQ32*
        Xb*Yt*(pow2(cbeta) + pow2(sbeta)))/((mD32 - mQ32)*pow2(sbeta)) + (pow2(
        cbeta)*pow2(Yt))/pow2(sbeta) + (2*(mD32 - 3*mQ32)*mQ32*Yt*(pow2(cbeta)
        + pow2(sbeta))*pow3(Xb))/(pow2(sbeta)*pow3(-mD32 + mQ32)) + (Xt2*pow2(
        sbeta)*pow2(-mD32 + mQ32 + Xb*Yb)*pow4(mQ3))/(pow2(cbeta)*pow2(mD32 -
        mQ32)*pow2(mQ32 - mU32)) + (Xb2*(2*pow2(cbeta)*pow2(sbeta)*pow4(mQ3) +
        pow4(cbeta)*(-2*mD32*pow2(Yt) + 3*mQ32*pow2(Yt) + pow4(mQ3)) + pow4(
        mQ3)*pow4(sbeta)))/(pow2(cbeta)*pow2(mD32 - mQ32)*pow2(sbeta)) + (pow2(
        cbeta)*pow2(Yt)*(-4*mD32*mQ32 + pow4(mD3) + 6*pow4(mQ3))*pow4(Xb))/(
        pow2(sbeta)*pow4(mD32 - mQ32)) + (2*mQ32*Xt*(mQ32*(-mD32 + mQ32)*Xb*(-
        mD32 + mQ32 + Xb*Yb)*pow2(sbeta) + pow2(cbeta)*(-3*mQ32*Yb*Yt*pow3(Xb)
        + (mQ32*(Xb - 3*Yt) - Xb*(Xb + Yb)*Yt)*pow4(mD3) + Xb*(Xb*Yb - 2*Xb*Yt
        - Yb*Yt)*pow4(mQ3) + mD32*(mQ32*Xb*(-(Xb*Yb) + 3*Xb*Yt + 2*Yb*Yt) + Yb*
        Yt*pow3(Xb) + (-2*Xb + 3*Yt)*pow4(mQ3)) + Yt*pow6(mD3) + (Xb - Yt)*
        pow6(mQ3))))/((mQ32 - mU32)*pow2(cbeta)*pow3(mD32 - mQ32))) + (pow2(
        cbeta)*pow2(-mD32 + mQ32 + Xb2)*pow2(Yt)*((-2*mQ32 + mU32)*pow2(-(mQ32*
        mU3) + pow3(mU3)) + pow4(mA)*(4*mQ32*mU32 + 5*pow4(mQ3) + 6*pow4(mU3))
        - 4*(mQ32 + mU32)*pow6(mA) - 2*mA2*(-(mU32*pow4(mQ3)) - 2*mQ32*pow4(
        mU3) + pow6(mQ3) + 2*pow6(mU3)) + pow8(mA)))/(deltaxyz(mA2,mU32,mQ32)*
        pow2(mD32 - mQ32)*pow2(sbeta))))/pow6(mQ3) + Xt2*((-12*Xb2*(3*mQ32*(-
        mD32 + mQ32)*(mD32 + mQ32 - 2*mU32)*dilog(1 - mD32/mQ32)*pow2(cbeta) +
        (mQ32 - mU32)*(3*mQ32*(-2*mD32 + mQ32 + mU32)*dilog(1 - mQ32/mU32)*
        pow2(cbeta) + (-mD32 + mQ32)*(-(mU32*pow2(cbeta)) + mQ32*(1 + pow2(
        cbeta)))) + 3*mQ32*dilog(1 - mD32/mU32)*pow2(cbeta)*pow2(mD32 - mU32)))
        /(mQ32*pow2(cbeta)*pow2(mD32 - mQ32)*pow2(mQ32 - mU32)) - (6*(4*mQ32*
        mU32 - 4*mQ32*mU32*dilog(1 - Mu2/mQ32) + 4*mQ32*mU32*dilog(1 - Mu2/
        mU32) - 6*mQ32*mU32*pow2(cbeta) + 6*mQ32*(mQ32 - mU32)*dilog(1 - mQ32/
        mU32)*pow2(cbeta) - 4*mQ32*mU32*pow2(sbeta) - 4*pow4(mQ3) + 5*pow2(
        cbeta)*pow4(mQ3) + 4*pow2(sbeta)*pow4(mQ3) + pow2(cbeta)*pow4(mU3)))/(
        mQ32*pow2(cbeta)*pow2(mQ32 - mU32)) - (6*pow4(Xb)*(-6*mQ32*(mQ32 -
        mU32)*dilog(1 - mD32/mU32)*pow2(mD32 - mU32) + 6*(mQ32 - mU32)*dilog(1
        - mD32/mQ32)*pow2(-(mD32*mQ3) + pow3(mQ3)) + 4*mQ32*mU32*pow4(mD3) - 6*
        mQ32*mU32*dilog(1 - mQ32/mU32)*pow4(mD3) - 32*mD32*mU32*pow4(mQ3) - 12*
        mD32*mU32*dilog(1 - mQ32/mU32)*pow4(mQ3) + 7*pow4(mD3)*pow4(mQ3) + 6*
        dilog(1 - mQ32/mU32)*pow4(mD3)*pow4(mQ3) + 10*mD32*mQ32*pow4(mU3) + 12*
        mD32*mQ32*dilog(1 - mQ32/mU32)*pow4(mU3) + pow4(mD3)*pow4(mU3) - 11*
        pow4(mQ3)*pow4(mU3) + 6*dilog(1 - mQ32/mU32)*pow4(mQ3)*pow4(mU3) - 2*
        mD32*pow6(mQ3) + 28*mU32*pow6(mQ3) - 6*mQ32*dilog(1 - mQ32/mU32)*pow6(
        mU3) - 5*pow8(mQ3)))/(pow2(-(mQ3*mU32) + pow3(mQ3))*pow4(mD32 - mQ32)));
}

double ThresholdCalculator::getDeltaLambdaYt4Yb2(int omitLogs) const
{
   using std::log;
   using himalaya::dilog;

   const double MR2 = pow2(p.scale);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mU32 = p.mu2(2,2);
   const double mQ3 = std::sqrt(mQ32);
   const double mU3 = std::sqrt(mU32);
   const double mD3 = std::sqrt(mD32);
   const double mA = p.MA;
   const double beta = std::atan(p.vu/p.vd);
   const double sbeta = std::sin(beta);
   const double cbeta = std::cos(beta);
   const double Xt = p.Au(2,2) - p.mu*p.vd/p.vu;
   const double Yt = p.Au(2,2) + p.mu*p.vu/p.vd;
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double Yb = p.Ad(2,2) + p.mu*p.vd/p.vu;
   const double Mu = p.mu;
   const double Mu2 = pow2(Mu);
   const double Xt2 = pow2(Xt);
   const double Xt4 = pow2(Xt2);
   const double mA2 = pow2(mA);
   const double Xb2 = pow2(Xb);
   const double lmQ3MR = omitLogs*log(mQ32 / MR2);
   const double lmUMR = omitLogs*log(Mu2 / MR2);
   const double eps = mQ3*0.01;

   if(std::abs(mU3 - mQ3) < eps && std::abs(mD3 - mQ3) < eps && std::abs(Mu - mQ3) < eps){
     const double M = std::sqrt(mQ3*mU3);
     return (2.675765602178715*pow2(sbeta)*(cbeta*M*sbeta*(-0.7474496265186756*Xb*Xt4
        + (8.96939551822411*Xb + 1.4948992530373513*Xt)*Xt2*pow2(M) - 
        8.96939551822411*Xb*pow4(M)) + pow4(cbeta)*(Yt*(1.4948992530373513*Xb*
        Xt2 + Xb2*(1.4948992530373513*Xt + 0.7474496265186756*Yt) + 
        0.7474496265186756*Xt2*Yt)*pow2(M) - 0.20235594112794902*Xb2*Xt2*pow2(
        Yt) + (1.2615562429128426*Xb2 - 7.0078102449377395*Xb*Xt + 
        1.2615562429128426*Xt2 - 7.0078102449377395*Xb*Yt - 7.0078102449377395*
        Xt*Yt - 3.5039051224688698*pow2(Yt))*pow4(M) - 10.511715367406609*pow6(
        M)) + pow2(sbeta)*(0.05195844587770295*Xt4*pow2(sbeta)*pow2(Yb) + Xt2*
        pow2(M)*(Xt2*(0.3737248132593378 - 0.560587219889007*pow2(sbeta)) - 
        0.342737744262778*Xt*Yb*pow2(sbeta) - 0.514106616394167*pow2(sbeta)*
        pow2(Yb)) + pow2(sbeta)*(5.746254002024897*Xt2 + 5.046224971651371*Xt*
        Yb + 1.2615562429128429*pow2(Yb))*pow4(M) + (8.969395518224111 + 1.*
        pow2(sbeta))*pow6(M)) + pow2(cbeta)*(Xt*pow2(M)*(Xb2*Xt*(
        1.1211744397780137 - 2.2423488795560274*pow2(sbeta)) + Xt2*(-
        1.1211744397780137*Xt - 0.342737744262778*Yb - 0.342737744262778*Yt)*
        pow2(sbeta) + Xb*(1.4948992530373513*Xt*Yb + 1.4948992530373513*Xt*Yt +
        2.9897985060747025*Yb*Yt)*pow2(sbeta)) + Xb*(0.3737248132593378*Xb*Xt -
        0.6125456657667098*Yb*Yt)*pow2(sbeta)*pow3(Xt) + (Xb*(-
        7.0078102449377395*Xt - 7.0078102449377395*Yb - 7.0078102449377395*Yt)*
        pow2(sbeta) + Xb2*(-2.2423488795560274 + 8.96939551822411*pow2(sbeta))
        + Xt*((5.0462249716513705*Yb + 5.0462249716513705*Yt)*pow2(sbeta) + Xt*
        (-2.2423488795560274 + 9.530922730763425*pow2(sbeta))))*pow4(M) + (-
        17.93879103644822 + 6.242348879556027*pow2(sbeta))*pow6(M))))/pow6(M);
   }

   return -4 + (48*Xb*std::abs(Mu)*((mD32 + Mu2)*(-mQ32 + Mu2)*dilog(1 - Mu2/mD32) + (
        mD32 - Mu2)*(mQ32 + Mu2)*dilog(1 - Mu2/mQ32)))/(cbeta*(mD32 - mQ32)*(
        mD32 - Mu2)*(mQ32 - Mu2)*sbeta) + (48*Xt*std::abs(Mu)*(dilog(1 - Mu2/mQ32) -
        dilog(1 - Mu2/mU32)))/(cbeta*(mQ32 - mU32)*sbeta) - 9/pow2(cbeta) - (6*
        mD32)/(mQ32*pow2(cbeta)) + (6*mD32)/((mD32 - Mu2)*pow2(cbeta)) + (12*
        Mu2)/(mQ32*pow2(cbeta)) + 8*pow2(Pi) + 12*dilog(1 - Mu2/mQ32)*(1/pow2(
        cbeta) - 2/pow2(sbeta)) - (24*dilog(1 - Mu2/mU32))/pow2(sbeta) + (3*
        pow2(sbeta))/pow2(cbeta) - (6*mA2*pow2(sbeta))/(mQ32*pow2(cbeta)) + (2*
        pow2(Pi)*pow2(sbeta))/pow2(cbeta) + 6*pow2(lmQ3MR)*(3 + (2*(-1 + pow2(
        cbeta)))/pow2(sbeta) + pow2(sbeta)/pow2(cbeta)) - (6*pow2(sbeta)*pow2(
        Yb))/(mQ32*pow2(cbeta)) + (24 + (6*pow2(sbeta))/pow2(cbeta))*pow2(log(
        mA2/MR2)) + (48*(cbeta*(-mQ32 + 2*Mu2 - mU32)*std::abs(Mu)*dilog(1 - Mu2/
        mQ32) + cbeta*(mQ32 - 2*Mu2 + mU32)*std::abs(Mu)*dilog(1 - Mu2/mU32) - 2*(
        mQ32 - mU32)*(2*cbeta*std::abs(Mu) + sbeta*(Yb + Yt)*pow2(cbeta) + Yb*pow3(
        sbeta)))*pow3(Xt))/(sbeta*pow2(cbeta)*pow3(mQ32 - mU32)) - (6*Xb2*(4*
        mD32*mQ32 + 4*mD32*mQ32*dilog(1 - Mu2/mD32) - 4*mD32*mQ32*dilog(1 -
        Mu2/mQ32) - 4*mD32*mQ32*pow2(cbeta) - 6*mD32*mQ32*pow2(sbeta) + 6*(mD32
        - mQ32)*mQ32*dilog(1 - mD32/mQ32)*pow2(sbeta) + pow2(sbeta)*pow4(mD3) -
        4*pow4(mQ3) + 4*pow2(cbeta)*pow4(mQ3) + 5*pow2(sbeta)*pow4(mQ3)))/(
        mQ32*pow2(mD32 - mQ32)*pow2(sbeta)) + 6*phixyz(mA2,mQ32,mQ32)*((2*mA2*
        Xb*Yb)/((mD32 - mQ32)*mQ32) - (2*mA2*(mD32 - mQ32 - Xb*Yb)*Yt*pow3(Xt))
        /((mD32 - mQ32)*pow2(-(mQ3*mU32) + pow3(mQ3))) + (mA2*(mA2 - 5*mQ32)*
        Xb2*pow2(cbeta))/((-mD32 + mQ32)*pow2(sbeta)*pow4(mQ3)) + (mA2*((mA2 -
        5*mQ32)*pow2(cbeta) - 2*mQ32*pow2(sbeta)))/(pow2(sbeta)*pow4(mQ3)) - (
        2*mA2*Xt*Yt*(-((mA2 - 5*mQ32)*(-mD32 + mQ32 + Xb2)*pow2(cbeta)) + mQ32*
        (-mD32 + mQ32 + Xb*Yb)*pow2(sbeta)))/((-mD32 + mQ32)*(mQ32 - mU32)*
        pow2(sbeta)*pow4(mQ3)) - (deltaxyz(mA2,mQ32,mQ32)*(2*Xt*pow2(mQ32 -
        mU32)*(-((mD32 - 2*mQ32)*Xb2*Yt*pow2(cbeta)) + Yt*pow2(cbeta)*pow2(mD32
        - mQ32) + (mD32 - mQ32)*mQ32*Xb*(pow2(cbeta) + pow2(sbeta))) - (mQ32 -
        mU32)*Xt2*(-2*(mD32 - mQ32)*mQ32*Xb*(Yt*pow2(cbeta) + (-Yb + Yt)*pow2(
        sbeta)) + (mD32 - 2*mQ32)*Xb2*pow2(cbeta)*pow2(Yt) - pow2(mD32 - mQ32)*
        (2*mQ32*pow2(sbeta) + pow2(cbeta)*pow2(Yt))) - (mD32 - 2*mQ32)*Xb2*
        pow2(cbeta)*pow3(mQ32 - mU32) + pow2(cbeta)*pow2(mD32 - mQ32)*pow3(mQ32
        - mU32) + 4*mQ32*(-mD32 + mQ32)*(-mD32 + mQ32 + Xb*Yb)*Yt*pow2(sbeta)*
        pow3(Xt)))/(pow2(mD32 - mQ32)*pow2(sbeta)*pow3(mQ32 - mU32)*pow4(mQ3))
        + (mA2*Xt2*(2*mQ32*(mQ32 - mU32)*Xb*Yb*pow2(sbeta) - (mA2 - 5*mQ32)*
        Xb2*pow2(cbeta)*pow2(Yt) - (mD32 - mQ32)*(-(mA2*pow2(cbeta)*pow2(Yt)) +
        mQ32*(-2*mU32*pow2(sbeta) + 5*pow2(cbeta)*pow2(Yt)) + 2*pow2(sbeta)*
        pow4(mQ3))))/((mD32 - mQ32)*pow2(mQ32 - mU32)*pow2(sbeta)*pow4(mQ3))) +
        6*pow2(log(mD32/mQ32))*((-4*(mD32 + Mu2)*Xb*std::abs(Mu))/(cbeta*(mD32 -
        mQ32)*(mD32 - Mu2)*sbeta) - (2*mD32*Xb2)/(pow2(mD32 - mQ32)*pow2(sbeta)
        ) + (-2*mD32*Mu2 + pow4(mD3) - pow4(Mu))/(pow2(cbeta)*pow2(mD32 - Mu2))
        ) + (12*dilog(1 - Mu2/mD32)*(-2*mD32*Mu2 + pow4(mD3) - pow4(Mu)))/(
        pow2(cbeta)*pow2(mD32 - Mu2)) + 6*lmUMR*(Mu2*(-((2/mQ32 + (2*mD32 +
        Mu2)/pow2(mD32 - Mu2))/pow2(cbeta)) + (4*(1/(mQ32 - Mu2) + 1/(-Mu2 +
        mU32)))/pow2(sbeta)) + (4*Mu2*Xb2)/((mD32 - mQ32)*(-mQ32 + Mu2)*pow2(
        sbeta)) + 4*Mu2*Xt2*(1/((mQ32 - mU32)*(Mu2 - mU32)*pow2(sbeta)) - 1/(
        pow2(cbeta)*(-(mQ32*mU32) + pow4(mQ3)))) + (2*Xt4*((4*cbeta*pow(Mu2,
        1.5)*(-mD32 + Mu2)*(mQ32 - mU32)*Xb)/((-mQ32 + Mu2)*sbeta) + (Mu2*(2*
        mQ32 + mU32)*pow4(mD3) - 2*mD32*(2*mQ32 + mU32)*pow4(Mu) + (Mu2*mU32 +
        mQ32*(2*Mu2 + mU32) - pow4(mQ3))*pow4(Mu))/mQ32))/(pow2(cbeta)*pow2(
        mD32 - Mu2)*pow3(mQ32 - mU32))) + 6*pow2(log(mU32/mQ32))*((-4*Xt*std::abs(
        Mu))/(cbeta*mQ32*sbeta - cbeta*mU32*sbeta) - 2/pow2(sbeta) + (Xt2*(2*(-
        Mu2 + mU32)*pow2(sbeta) + pow2(cbeta)*(3*(-2*mU32 + Xb2)*pow2(sbeta) +
        mQ32*(2 + 6*pow2(sbeta)))))/(pow2(cbeta)*pow2(mQ32 - mU32)*pow2(sbeta))
        + (4*(mQ32 - 2*Mu2 + mU32)*std::abs(Mu)*pow3(Xt))/(cbeta*sbeta*pow3(mQ32 -
        mU32)) + (Xt4*((Mu2 - mU32)*(2*mQ32 - 3*Mu2 + mU32) + 3*pow2(cbeta)*(
        mD32*Xb2 - 2*mU32*Xb2 + mQ32*(-2*mU32 + Xb2) + pow4(mQ3) + pow4(mU3))))
        /(pow2(cbeta)*pow4(mQ32 - mU32))) + 3*phixyz(mA2,mQ32,mD32)*((2*Xt2*Yb*
        ((-2*(mA2 + mD32 - mQ32)*(mQ32 - mU32)*Xb)/(mD32 - mQ32) + ((mA2 + mD32
        - mU32)*Yb*pow2(sbeta))/pow2(cbeta)))/(mD32*pow2(mQ32 - mU32)) - (4*(
        mA2 + mD32 - mQ32)*Xt*Yb*(Xb*Yt*pow2(cbeta) + (mD32 - mQ32)*(pow2(
        cbeta) + pow2(sbeta))))/(mD32*(mD32 - mQ32)*(mQ32 - mU32)*pow2(cbeta))
        + (pow2(sbeta)*pow2(Yb))/(mD32*pow2(cbeta)) - (pow2(mA2 + mD32 - mQ32)*
        pow2(sbeta)*pow2(mQ32 - mU32 + Xt2)*pow2(Yb))/(mD32*deltaxyz(mA2,mQ32,
        mD32)*pow2(cbeta)*pow2(mQ32 - mU32)) + ((4*mA2 + 4*mD32 - 3*mQ32 -
        mU32)*Xt4*pow2(sbeta)*pow2(Yb))/(mD32*pow2(cbeta)*pow3(mQ32 - mU32)) -
        (4*(mA2 + mD32 - mQ32)*Yb*(Xb*Yt*pow2(cbeta) + (mD32 - mQ32)*(pow2(
        cbeta) + pow2(sbeta)))*pow3(Xt))/(mD32*(mD32 - mQ32)*pow2(cbeta)*pow2(
        mQ32 - mU32)) - (4*(mA2 + mD32 - mQ32)*Xb*Yb)/(-(mD32*mQ32) + pow4(mD3)
        ) + (2*deltaxyz(mA2,mQ32,mD32)*(2*Xb*Xt*pow2(cbeta)*(Xb*Yt*pow2(cbeta)
        + (mD32 - mQ32)*(pow2(cbeta) + pow2(sbeta)))*pow3(mQ32 - mU32) - 4*(
        mD32 - mQ32)*(mQ32 - mU32)*Yb*pow2(sbeta)*(Xb*Yt*pow2(cbeta) + (mD32 -
        mQ32)*(pow2(cbeta) + pow2(sbeta)))*pow3(Xt) + Xt2*pow2(mQ32 - mU32)*(2*
        (mD32 - mQ32)*Xb*pow2(cbeta)*(Yt*pow2(cbeta) + (-Yb + Yt)*pow2(sbeta))
        + pow2(mD32 - mQ32)*pow2(pow2(cbeta) + pow2(sbeta)) + Xb2*pow2(Yt)*
        pow4(cbeta)) + Xb2*pow4(cbeta)*pow4(mQ32 - mU32) + 3*Xt4*pow2(mD32 -
        mQ32)*pow2(Yb)*pow4(sbeta)))/(mD32*pow2(cbeta)*pow2(mD32 - mQ32)*pow2(
        sbeta)*pow4(mQ32 - mU32))) + (12*Xt2*((Xb2*(3*mQ32*(-mD32 + mQ32)*(mD32
        + mQ32 - 2*mU32)*dilog(1 - mD32/mQ32)*pow2(sbeta) + 3*mQ32*dilog(1 -
        mD32/mU32)*pow2(mD32 - mU32)*pow2(sbeta) - (mQ32 - mU32)*(-3*mQ32*(-2*
        mD32 + mQ32 + mU32)*dilog(1 - mQ32/mU32)*pow2(sbeta) + (mD32 - mQ32)*(
        mD32*pow2(sbeta) - mQ32*(1 + pow2(sbeta))))))/pow2(mD32 - mQ32) + (2*
        mQ32*mU32*pow2(cbeta) - mD32*mQ32*pow2(sbeta) + 2*mQ32*Mu2*pow2(sbeta)
        + mD32*mU32*pow2(sbeta) + 3*mQ32*mU32*pow2(sbeta) - 2*Mu2*mU32*pow2(
        sbeta) - 2*mQ32*Mu2*dilog(1 - Mu2/mU32)*pow2(sbeta) + 2*mQ32*mU32*
        dilog(1 - Mu2/mU32)*pow2(sbeta) - 4*mQ32*mU32*pow2(cbeta)*pow2(sbeta) +
        6*mQ32*(mQ32 - mU32)*dilog(1 - mQ32/mU32)*pow2(cbeta)*pow2(sbeta) - 2*
        mQ32*mU32*pow4(cbeta) - 2*pow2(cbeta)*pow4(mQ3) + 2*dilog(1 - Mu2/mU32)
        *pow2(cbeta)*pow4(mQ3) - 3*pow2(sbeta)*pow4(mQ3) + 4*pow2(cbeta)*pow2(
        sbeta)*pow4(mQ3) + 2*pow4(cbeta)*pow4(mQ3) - 2*dilog(1 - Mu2/mQ32)*(
        mQ32*(-Mu2 + mU32)*pow2(sbeta) + pow2(cbeta)*pow4(mQ3)) - mA2*mQ32*
        pow4(sbeta) + mA2*mU32*pow4(sbeta) - 2*mQ32*mU32*pow4(sbeta) - mQ32*
        pow2(Yb)*pow4(sbeta) + mU32*pow2(Yb)*pow4(sbeta) + 2*pow4(mQ3)*pow4(
        sbeta))/pow2(cbeta)))/(mQ32*pow2(mQ32 - mU32)*pow2(sbeta)) + (3*log(
        mA2/MR2)*(-4*lmQ3MR*(4*pow2(cbeta) + pow2(sbeta)) + (2*Xb2*pow2(Yt)*
        pow4(cbeta))/((mD32 - mQ32)*mQ32*pow2(sbeta)) + (2*(-mD32 + mQ32 + Xb2)
        *(mQ32 - mU32 - Xt2)*pow2(mA2 + mQ32 - mU32)*pow2(Yt)*pow4(cbeta))/(
        mQ32*(-mD32 + mQ32)*(mQ32 - mU32)*deltaxyz(mA2,mU32,mQ32)*pow2(sbeta))
        + (pow2(sbeta)*pow2(mQ32 - mU32 + Xt2)*pow2(Yb)*(-pow2(mD32 - mQ32) +
        pow4(mA)))/(deltaxyz(mA2,mQ32,mD32)*pow2(-(mQ3*mU32) + pow3(mQ3))) - (
        Xt4*(pow2(sbeta)*(2*mA2*(5*mQ32 + mU32) + (mQ32 - mU32)*(4*mQ32 - pow2(
        Yb))) + 16*pow2(cbeta)*(-(mQ32*mU32) + pow4(mQ3))))/(mQ32*pow3(mQ32 -
        mU32)) + ((8*mQ32 - 2*pow2(Yt))*pow4(cbeta) + (2*mA2 + 2*mQ32 + pow2(
        Yb))*pow4(sbeta))/(mQ32*pow2(sbeta)) - (2*Xt2*(Xb2*pow2(Yt)*pow4(cbeta)
        - (mD32 - mQ32)*(pow2(Yt)*pow4(cbeta) + 2*mA2*pow4(sbeta) + pow2(Yb)*
        pow4(sbeta))))/((mD32 - mQ32)*mQ32*(mQ32 - mU32)*pow2(sbeta))))/pow2(
        cbeta) - (6*phixyz(mA2,mU32,mD32)*(-2*(mD32 - mQ32)*(mA2 + mD32 - mU32)
        *Xt*Yb*pow2(sbeta)*(Xb*Yt*pow2(cbeta) + (mD32 - mQ32)*(pow2(cbeta) +
        pow2(sbeta)))*pow3(mQ32 - mU32) + 2*(mD32 - mQ32)*(mA2 + mD32 - mU32)*
        Yb*pow2(mQ32 - mU32)*pow2(sbeta)*(Xb*Yt*pow2(cbeta) + (mD32 - mQ32)*(
        pow2(cbeta) + pow2(sbeta)))*pow3(Xt) + 2*(mD32 - mQ32)*(mA2 + mD32 -
        mU32)*Xb*Yt*pow2(cbeta)*(pow2(cbeta) + pow2(sbeta))*pow4(mQ32 - mU32) +
        (mA2 + mD32 - mU32)*pow2(mD32 - mQ32)*pow2(pow2(cbeta) + pow2(sbeta))*
        pow4(mQ32 - mU32) + (mA2 + mD32 - mU32)*Xb2*pow2(Yt)*pow4(cbeta)*pow4(
        mQ32 - mU32) - (mA2 + mD32 - mU32)*(mQ32 - mU32)*Xt4*pow2(mD32 - mQ32)*
        pow2(Yb)*pow4(sbeta) + Xt*deltaxyz(mA2,mU32,mD32)*(-4*(mD32 - mQ32)*(
        mQ32 - mU32)*Xt2*Yb*pow2(sbeta)*(Xb*Yt*pow2(cbeta) + (mD32 - mQ32)*(
        pow2(cbeta) + pow2(sbeta))) + 2*Xb*pow2(cbeta)*(Xb*Yt*pow2(cbeta) + (
        mD32 - mQ32)*(pow2(cbeta) + pow2(sbeta)))*pow3(mQ32 - mU32) + Xt*pow2(
        mQ32 - mU32)*(2*(mD32 - mQ32)*Xb*pow2(cbeta)*(Yt*pow2(cbeta) + (-Yb +
        Yt)*pow2(sbeta)) + pow2(mD32 - mQ32)*pow2(pow2(cbeta) + pow2(sbeta)) +
        Xb2*pow2(Yt)*pow4(cbeta)) + 3*pow2(mD32 - mQ32)*pow2(Yb)*pow3(Xt)*pow4(
        sbeta)) + (mA2 + mD32 - mU32)*Xt2*pow2(mQ32 - mU32)*(-2*(mD32 - mQ32)*(
        mQ32 - mU32)*Xb*Yt*pow2(cbeta)*(pow2(cbeta) + pow2(sbeta)) - (mQ32 -
        mU32)*Xb2*pow2(Yt)*pow4(cbeta) + pow2(mD32 - mQ32)*(2*(-mQ32 + mU32)*
        pow2(cbeta)*pow2(sbeta) + (-mQ32 + mU32)*pow4(cbeta) + (-mQ32 + mU32 +
        pow2(Yb))*pow4(sbeta)))))/(mD32*pow2(cbeta)*pow2(mD32 - mQ32)*pow2(
        sbeta)*pow4(mQ32 - mU32)) + 6*phixyz(mA2,mU32,mQ32)*((2*(mA2 + mQ32 -
        mU32)*Xb*Yt*(pow2(cbeta) + pow2(sbeta)))/((mD32 - mQ32)*mQ32*pow2(
        sbeta)) + (2*(mA2 + mQ32 - mU32)*(-mD32 + mQ32 + Xb*Yb)*Yt*pow3(Xt))/((
        mD32 - mQ32)*pow2(-(mQ3*mU32) + pow3(mQ3))) + ((mA2 - mU32)*pow2(cbeta)
        *pow2(Yt))/(pow2(sbeta)*pow4(mQ3)) + (Xt*deltaxyz(mA2,mU32,mQ32)*(4*
        mQ32*(-mD32 + mQ32)*Xt2*(-mD32 + mQ32 + Xb*Yb)*Yt*pow2(sbeta) + 2*pow2(
        mQ32 - mU32)*(-((mD32 - 2*mQ32)*Xb2*Yt*pow2(cbeta)) + Yt*pow2(cbeta)*
        pow2(mD32 - mQ32) + (mD32 - mQ32)*mQ32*Xb*(pow2(cbeta) + pow2(sbeta)))
        - (mQ32 - mU32)*Xt*(-2*(mD32 - mQ32)*mQ32*Xb*(Yt*pow2(cbeta) + (-Yb +
        Yt)*pow2(sbeta)) + (mD32 - 2*mQ32)*Xb2*pow2(cbeta)*pow2(Yt) - pow2(mD32
        - mQ32)*(2*mQ32*pow2(sbeta) + pow2(cbeta)*pow2(Yt)))))/(pow2(mD32 -
        mQ32)*pow2(sbeta)*pow3(mQ32 - mU32)*pow4(mQ3)) + (Xb2*pow2(cbeta)*pow2(
        Yt)*(-(mA2*(mD32 - 2*mQ32)) + mD32*mU32 - 2*mQ32*mU32 + pow4(mQ3)))/(
        pow2(mD32 - mQ32)*pow2(sbeta)*pow4(mQ3)) - (Xt2*Yt*(2*(mD32 - mQ32)*
        mQ32*(mQ32 - mU32)*(mA2 + mQ32 - mU32)*Xb*(pow2(cbeta) + pow2(sbeta)) +
        Yt*pow2(cbeta)*pow2(mD32 - mQ32)*(-(mA2*(2*mQ32 + 3*mU32)) + 2*pow2(
        mQ32 - mU32) + pow4(mA)) - Xb2*Yt*pow2(cbeta)*((2*mD32 - 3*mQ32)*pow2(
        mQ32 - mU32) + (mD32 - mQ32)*pow4(mA) + mA2*(4*mQ32*mU32 - mD32*(2*mQ32
        + 3*mU32) + pow4(mQ3)))))/(pow2(mD32 - mQ32)*pow2(mQ32 - mU32)*pow2(
        sbeta)*pow4(mQ3)) + (2*Xt*Yt*(mQ32*(mA2 + mQ32 - mU32)*(-mD32 + mQ32 +
        Xb*Yb)*pow2(sbeta) - (-mD32 + mQ32 + Xb2)*pow2(cbeta)*(-3*mQ32*mU32 -
        mA2*(3*mQ32 + 2*mU32) + pow4(mA) + 2*pow4(mQ3) + pow4(mU3))))/((-mD32 +
        mQ32)*(mQ32 - mU32)*pow2(sbeta)*pow4(mQ3)) + ((-mD32 + mQ32 + Xb2)*(
        mQ32 - mU32 - Xt2)*pow2(cbeta)*pow2(Yt)*(pow2(-(mQ32*mU3) + pow3(mU3))
        + 3*mU32*pow4(mA) + mA2*(2*mQ32*mU32 + pow4(mQ3) - 3*pow4(mU3)) - pow6(
        mA)))/((-mD32 + mQ32)*(mQ32 - mU32)*deltaxyz(mA2,mU32,mQ32)*pow2(sbeta)
        *pow4(mQ3))) + 3*lmQ3MR*(4 + (8*Mu2*(1/(-mQ32 + Mu2) + 1/(Mu2 - mU32)))
        /pow2(sbeta) - (2*(-mD32 + mQ32 + Xb2)*(mQ32 - mU32 - Xt2)*pow2(cbeta)*
        pow2(mA2 + mQ32 - mU32)*pow2(Yt))/(mQ32*(-mD32 + mQ32)*(mQ32 - mU32)*
        deltaxyz(mA2,mU32,mQ32)*pow2(sbeta)) + (2*pow2(cbeta)*(-4 + pow2(Yt)/
        mQ32))/pow2(sbeta) + (16*(2*cbeta*std::abs(Mu) + sbeta*(Yb + Yt)*pow2(cbeta)
        + Yb*pow3(sbeta))*pow3(Xt))/(sbeta*pow2(cbeta)*pow2(mQ32 - mU32)) + (2*
        Xt2*((Xb2*(2*mD32*pow2(sbeta) - 2*mQ32*pow2(sbeta) + pow2(cbeta)*pow2(
        Yt)))/(mD32 - mQ32) + (2*mQ32*pow2(cbeta)*(mU32*(1 - 2*pow2(sbeta)) +
        Mu2*(1 + 2*pow2(sbeta))) - (Mu2 - mU32)*pow2(sbeta)*(2*mD32 - 2*mQ32*(-
        2 + pow2(sbeta)) + pow2(sbeta)*pow2(Yb)) + (Mu2 - mU32)*(2*mQ32 + pow2(
        Yt))*pow4(cbeta))/((-Mu2 + mU32)*pow2(cbeta))))/(mQ32*(mQ32 - mU32)*
        pow2(sbeta)) + (pow2(sbeta)*pow2(mQ32 - mU32 + Xt2)*pow2(Yb)*(pow2(mD32
        - mQ32) - pow4(mA)))/(mQ32*deltaxyz(mA2,mQ32,mD32)*pow2(cbeta)*pow2(
        mQ32 - mU32)) + (4 + mD32*(2/mQ32 + (8*Mu2)/pow2(mD32 - Mu2)) + pow2(
        sbeta)*(-2 + pow2(Yb)/mQ32) - (2*pow4(mD3))/pow2(mD32 - Mu2))/pow2(
        cbeta) + (2*Xb2*(Mu2*(mD32*pow2(sbeta) - pow2(cbeta)*pow2(Yt)) - mQ32*(
        mD32*pow2(sbeta) + Mu2*(2 + 2*pow2(cbeta) + 3*pow2(sbeta)) - pow2(
        cbeta)*pow2(Yt)) + (-2 + 2*pow2(cbeta) + 3*pow2(sbeta))*pow4(mQ3)))/(
        mQ32*(-mD32 + mQ32)*(mQ32 - Mu2)*pow2(sbeta)) + (Xt4*((-16*pow(Mu2,
        1.5)*(mQ32 - mU32)*Xb)/(cbeta*(-mD32 + Mu2)*(-mQ32 + Mu2)*sbeta) - (2*(
        5*mQ32 + mU32)*Xb2)/mQ32 - (pow4(mD3)*(mU32*(-4*Mu2 + pow2(sbeta)*pow2(
        Yb)) + mQ32*(-32*Mu2 + mU32*(6 - 16*pow2(cbeta)) + 11*pow2(sbeta)*pow2(
        Yb)) + 2*(3 + 8*pow2(cbeta))*pow4(mQ3)) - 2*mD32*Mu2*(mU32*(-Mu2 +
        pow2(sbeta)*pow2(Yb)) + mQ32*(-17*Mu2 + mU32*(6 - 16*pow2(cbeta)) + 11*
        pow2(sbeta)*pow2(Yb)) + 2*(3 + 8*pow2(cbeta))*pow4(mQ3)) + (mU32*pow2(
        sbeta)*pow2(Yb) + mQ32*(-12*Mu2 - 2*mU32*(-5 + 8*pow2(cbeta)) + 11*
        pow2(sbeta)*pow2(Yb)) + 2*(1 + 8*pow2(cbeta))*pow4(mQ3))*pow4(Mu) + 2*(
        5*mQ32 + mU32)*pow6(mD3))/(mQ32*pow2(cbeta)*pow2(mD32 - Mu2))))/pow3(
        mQ32 - mU32)) + log(mU32/mQ32)*((9 + (6*mD32)/(-mD32 + Mu2) - 9*pow2(
        sbeta))/pow2(cbeta) + (24*Mu2)/((Mu2 - mU32)*pow2(sbeta)) - (6*pow2(
        cbeta)*(2*mQ32 - pow2(Yt)))/(mQ32*pow2(sbeta)) + (6*Xb2*(2*mQ32 + pow2(
        cbeta)*pow2(Yt)))/(mQ32*(-mD32 + mQ32)*pow2(sbeta)) + (24*Xt*(2*cbeta*
        std::abs(Mu) + sbeta*(Yb + Yt)*pow2(cbeta) + Yb*pow3(sbeta)))/((mQ32 - mU32)
        *sbeta*pow2(cbeta)) - (24*(mQ32 + 3*mU32)*(2*cbeta*std::abs(Mu) + sbeta*(Yb
        + Yt)*pow2(cbeta) + Yb*pow3(sbeta))*pow3(Xt))/(sbeta*pow2(cbeta)*pow3(
        mQ32 - mU32)) - (6*(-mD32 + mQ32 + Xb2)*(mQ32 - mU32 - Xt2)*pow2(cbeta)
        *pow2(Yt)*(-2*mA2*mU32 + pow4(mA) - pow4(mQ3) + pow4(mU3)))/(mQ32*(-
        mD32 + mQ32)*(mQ32 - mU32)*deltaxyz(mA2,mU32,mQ32)*pow2(sbeta)) + (6*
        Xt2*((mD32*(6*Mu2 - 5*mU32 + mQ32*(-1 + pow2(sbeta)) - 2*mA2*pow2(
        sbeta) + 3*mU32*pow2(sbeta) - 2*pow2(sbeta)*pow2(Yb)) + Mu2*(-4*Mu2 +
        3*mU32 - mQ32*(-3 + pow2(sbeta)) + 2*mA2*pow2(sbeta) - 3*mU32*pow2(
        sbeta) + 2*pow2(sbeta)*pow2(Yb)) - 2*pow4(mD3))/((mD32 - Mu2)*pow2(
        cbeta)) + (pow2(cbeta)*(mQ32*(2*mU32 - pow2(Yt)) + mU32*pow2(Yt) + 2*
        pow4(mQ3)))/(mQ32*pow2(sbeta)) - (Xb2*(-(mU32*pow2(cbeta)*pow2(Yt)) +
        mQ32*(-2*mD32*pow2(sbeta) + 2*mU32*pow2(sbeta) + pow2(cbeta)*pow2(Yt))
        + 2*pow4(mQ3)))/(mQ32*(-mD32 + mQ32)*pow2(sbeta)) + (2*(7*Mu2*mU32*
        pow2(sbeta) + mQ32*(3*mU32*pow2(sbeta) - Mu2*(2 + 3*pow2(sbeta))) + (2
        - 7*pow2(sbeta))*pow4(mU3)))/((Mu2 - mU32)*pow2(sbeta))))/pow2(mQ32 -
        mU32) + (3*Xt4*((4*(mQ32 + 2*mU32)*pow4(mD3) + mD32*(4*mQ32*(-3*Mu2 +
        mU32*(3 + 6*pow2(cbeta)) + pow2(sbeta)*(mA2 + pow2(Yb))) + mU32*(-36*
        Mu2 - mU32*(-11 + 36*pow2(cbeta) + pow2(sbeta)) + 4*pow2(sbeta)*(2*mA2
        + 5*pow2(Yb))) + (1 + 12*pow2(cbeta) + pow2(sbeta))*pow4(mQ3)) + Mu2*(
        4*mQ32*(2*Mu2 - 3*mU32*(1 + 2*pow2(cbeta)) - pow2(sbeta)*(mA2 + pow2(
        Yb))) + mU32*(28*Mu2 + mU32*(-9 + 36*pow2(cbeta) + pow2(sbeta)) - 4*
        pow2(sbeta)*(2*mA2 + 5*pow2(Yb))) - (3 + 12*pow2(cbeta) + pow2(sbeta))*
        pow4(mQ3)))/((mD32 - Mu2)*pow2(cbeta)) + (4*Xb2*(-7*mQ32*mU32 + mD32*(
        mQ32 + 5*mU32) + pow4(mU3)))/(mD32 - mQ32)))/pow4(mQ32 - mU32) + (6*
        lmUMR*(16*cbeta*(mQ32 - Mu2)*pow(Mu2,1.5)*(mQ32 - mU32)*pow2(mD32 -
        Mu2)*pow3(Xt) + 2*Xt2*pow3(mQ32 - mU32)*(4*cbeta*pow(Mu2,1.5)*(-mD32
        + Mu2)*Xb + (mQ32 - Mu2)*sbeta*pow4(Mu)) + Xt4*(4*cbeta*(mD32 - Mu2)*
        pow(Mu2,1.5)*(mQ32 - mU32)*(mQ32 + mU32)*Xb + (mQ32 - Mu2)*sbeta*
        pow4(Mu)*(-12*mD32*Mu2 + 6*pow4(mD3) - pow4(mQ3) + 6*pow4(Mu) + pow4(
        mU3))) + 4*cbeta*(mD32 - Mu2)*pow(Mu2,1.5)*Xb*pow4(mQ32 - mU32) + (-
        mQ32 + Mu2)*sbeta*pow4(Mu)*pow4(mQ32 - mU32)))/((mQ32 - Mu2)*sbeta*
        pow2(cbeta)*pow2(mD32 - Mu2)*pow4(mQ32 - mU32)) + 6*log(mA2/MR2)*(-2 +
        pow2(cbeta)/pow2(sbeta) + (2*((-2*mD32*Yb + mQ32*Yb + mU32*Yb - mQ32*Yt
        + mU32*Yt - 2*Xb*Yb*Yt + 2*mA2*(Yb + Yt))*pow2(cbeta) + (2*mA2 - 2*mD32
        + mQ32 + mU32)*Yb*pow2(sbeta))*pow3(Xt))/(pow2(cbeta)*pow3(mQ32 - mU32)
        ) - (Xt4*(4*pow2(cbeta)*(pow4(mQ3) - pow4(mU3)) + pow2(sbeta)*(-3*mD32*
        pow2(Yb) + mQ32*pow2(Yb) + 2*mU32*pow2(Yb) + mA2*(2*mQ32 + 4*mU32 + 3*
        pow2(Yb)) + pow4(mQ3) - pow4(mU3))))/(pow2(cbeta)*pow4(mQ32 - mU32)) +
        (2*Xt*((Xb - Yb - Yt)*pow2(cbeta)*pow2(sbeta) + (Xb + Yt)*pow4(cbeta) -
        Yb*pow4(sbeta)))/((mQ32 - mU32)*pow2(cbeta)*pow2(sbeta)) - (Xt2*(-2*(
        mD32 + 2*mQ32 - 3*mU32 - Xb*Yb + Xb*Yt)*pow2(cbeta)*pow2(sbeta) + (mA2
        - mD32 + mQ32 - 2*Xb*Yt - pow2(Yt))*pow4(cbeta) - (mA2 + mD32 + mQ32 -
        2*mU32 + pow2(Yb))*pow4(sbeta)))/(pow2(cbeta)*pow2(mQ32 - mU32)*pow2(
        sbeta))) + 6*lmQ3MR*(2 - (4*pow(Mu2,1.5)*Xb)/(cbeta*(-mD32 + Mu2)*(-
        mQ32 + Mu2)*sbeta) - 2/pow2(sbeta) + pow2(cbeta)/pow2(sbeta) + (2*(4*
        cbeta*(mQ32 - 2*Mu2 + mU32)*std::abs(Mu) + sbeta*(2*mD32*Yb + mQ32*Yb +
        mU32*Yb + 3*mQ32*Yt + mU32*Yt + 2*Xb*Yb*Yt - 2*mA2*(Yb + Yt))*pow2(
        cbeta) + (-2*mA2 + 2*mD32 + mQ32 + mU32)*Yb*pow3(sbeta))*pow3(Xt))/(
        sbeta*pow2(cbeta)*pow3(mQ32 - mU32)) + ((2*mD32*Mu2)/pow2(mD32 - Mu2) +
        pow2(sbeta) - pow4(mD3)/pow2(mD32 - Mu2))/pow2(cbeta) - (2*Xt*(4*cbeta*
        sbeta*std::abs(Mu) + (Xb + Yb + Yt)*pow2(cbeta)*pow2(sbeta) + (Xb + Yt)*
        pow4(cbeta) + Yb*pow4(sbeta)))/((mQ32 - mU32)*pow2(cbeta)*pow2(sbeta))
        + (Xt2*(-2*mD32 - 6*mU32 + 2*Xb2 + mQ32*(4 + 2/pow2(sbeta)) - (pow2(
        cbeta)*(-mA2 + mD32 + mQ32 + pow2(Yt)))/pow2(sbeta) + (2*Xb*(4*pow(
        Mu2,1.5)*mU32*sbeta + cbeta*Mu2*(-mD32 + Mu2)*(Yt*pow2(cbeta) + (-Yb +
        Yt)*pow2(sbeta)) - mQ32*(4*pow(Mu2,1.5)*sbeta + cbeta*(mD32 - Mu2)*(
        Yb - Yt)*pow2(sbeta) + (-mD32 + Mu2)*Yt*pow3(cbeta))))/(cbeta*(mD32 -
        Mu2)*(-mQ32 + Mu2)*pow2(sbeta)) + (-(mD32*Mu2*(4*mU32 + Mu2*(-10 +
        pow2(sbeta)) - 2*mQ32*(-2 + pow2(sbeta)) + 2*pow2(sbeta)*(mA2 + pow2(
        Yb)))) + (2*mU32 + 2*Mu2*(-4 + pow2(sbeta)) - mQ32*(-2 + pow2(sbeta)) +
        mA2*pow2(sbeta) + pow2(sbeta)*pow2(Yb))*pow4(mD3) + (-4*Mu2 + 4*mU32 +
        pow2(sbeta)*(mA2 - mQ32 + pow2(Yb)))*pow4(Mu) - (-2 + pow2(sbeta))*
        pow6(mD3))/(pow2(cbeta)*pow2(mD32 - Mu2))))/pow2(mQ32 - mU32) + (Xt4*((
        -4*pow(Mu2,1.5)*(mQ32 - mU32)*(mQ32 + mU32)*Xb)/(cbeta*(-mD32 + Mu2)*
        (-mQ32 + Mu2)*sbeta) - 2*(mQ32 + 2*mU32)*Xb2 - (pow4(Mu)*(-8*Mu2*mU32 -
        3*mA2*pow2(sbeta)*pow2(Yb) + 2*mU32*pow2(sbeta)*pow2(Yb) + mQ32*(-4*Mu2
        + 4*mU32 + pow2(sbeta)*pow2(Yb)) + 6*pow4(Mu) + 4*pow2(cbeta)*(pow4(
        mQ3) - pow4(mU3)) + 2*pow4(mU3)) - mD32*Mu2*(-6*mA2*pow2(sbeta)*pow2(
        Yb) + 4*mU32*pow2(sbeta)*pow2(Yb) + mQ32*(-10*Mu2 + 8*mU32 + 2*pow2(
        sbeta)*pow2(Yb)) - Mu2*(20*mU32 + 3*pow2(sbeta)*pow2(Yb)) + (2 + 8*
        pow2(cbeta))*pow4(mQ3) + 12*pow4(Mu) + 2*pow4(mU3) - 8*pow2(cbeta)*
        pow4(mU3)) + pow4(mD3)*(-3*mA2*pow2(sbeta)*pow2(Yb) + 2*mU32*pow2(
        sbeta)*pow2(Yb) + mQ32*(-8*Mu2 + 4*mU32 + pow2(sbeta)*pow2(Yb)) - 2*
        Mu2*(8*mU32 + 3*pow2(sbeta)*pow2(Yb)) + (1 + 4*pow2(cbeta))*pow4(mQ3) +
        6*pow4(Mu) + pow4(mU3) - 4*pow2(cbeta)*pow4(mU3)) + (2*mQ32 + 4*mU32 +
        3*pow2(sbeta)*pow2(Yb))*pow6(mD3))/(pow2(cbeta)*pow2(mD32 - Mu2))))/
        pow4(mQ32 - mU32))) + log(mD32/mQ32)*((6*Xt2*((2*mD32*Xb2*(mD32*pow2(
        sbeta) - mQ32*(1 + pow2(sbeta))))/(pow2(mD32 - mQ32)*pow2(sbeta)) + (2*
        mD32 + pow2(sbeta)*pow2(Yb))/pow2(cbeta)))/(mQ32*(mQ32 - mU32)) + (3*(
        mD32*(2/mQ32 + (4*Mu2)/pow2(mD32 - Mu2)) + (pow2(sbeta)*pow2(Yb))/mQ32
        + (2*pow4(mD3))/pow2(mD32 - Mu2)))/pow2(cbeta) + (6*Xb2*(mD32*mQ32*(4 -
        4*pow2(cbeta) - 5*pow2(sbeta)) + pow2(sbeta)*pow4(mD3)))/(mQ32*pow2(
        mD32 - mQ32)*pow2(sbeta)) - (3*pow2(sbeta)*pow2(mQ32 - mU32 + Xt2)*
        pow2(Yb)*(-2*mA2*mQ32 + pow4(mA) - pow4(mD3) + pow4(mQ3)))/(mQ32*
        deltaxyz(mA2,mQ32,mD32)*pow2(cbeta)*pow2(mQ32 - mU32)) + (12*lmUMR*((4*
        cbeta*(mD32 - Mu2)*pow(Mu2,1.5)*Xb)/((mD32 - mQ32)*sbeta) + pow4(Mu))
        )/(pow2(cbeta)*pow2(mD32 - Mu2)) + 6*lmQ3MR*(2 + pow2(cbeta)/pow2(
        sbeta) + (2*Xb*Xt*((mD32 - mQ32 + Xb*Yt)*pow2(cbeta) + (mD32 - mQ32)*
        pow2(sbeta)))/(pow2(mD32 - mQ32)*pow2(sbeta)) + (Xb2*(2*mD32*(-1 +
        pow2(sbeta)) + pow2(cbeta)*(-mA2 + mD32 + mQ32 + pow2(Yt))))/(pow2(mD32
        - mQ32)*pow2(sbeta)) + (2*Xb*(-4*pow(Mu2,1.5)*sbeta + cbeta*(mD32 -
        Mu2)*(Yb + Yt)*pow2(sbeta) + (mD32 - Mu2)*Yt*pow3(cbeta)))/(cbeta*(mD32
        - mQ32)*(mD32 - Mu2)*pow2(sbeta)) + (pow2(sbeta) - (2*pow4(Mu))/pow2(
        mD32 - Mu2))/pow2(cbeta)) - (6*log(mA2/MR2)*(2*(mD32 - mQ32)*(mD32 -
        mQ32 + Xb*(Xt + Yb + Yt))*pow2(cbeta)*pow2(sbeta) + pow4(cbeta)*(mQ32*
        Xb*(Xb - 2*(Xt + Yt)) + Xb2*(-mA2 + Yt*(2*Xt + Yt)) - mD32*(2*mQ32 +
        Xb*(Xb - 2*(Xt + Yt))) + pow4(mD3) + pow4(mQ3)) + pow2(mD32 - mQ32)*
        pow4(sbeta)))/(pow2(cbeta)*pow2(mD32 - mQ32)*pow2(sbeta)) + Xt4*((48*
        mD32*Xb*std::abs(Mu))/(cbeta*(mD32 - mQ32)*(mD32 - Mu2)*sbeta*pow2(mQ32 -
        mU32)) + (6*mD32*Xb2)/((mD32 - mQ32)*pow2(-(mQ3*mU32) + pow3(mQ3))) - (
        3*(2*mD32*Mu2*(mU32*(Mu2 - pow2(sbeta)*pow2(Yb)) + mQ32*(5*Mu2 + 4*mU32
        + pow2(sbeta)*pow2(Yb)) - 4*pow4(mQ3)) + pow4(mD3)*(mU32*(-4*Mu2 +
        pow2(sbeta)*pow2(Yb)) - mQ32*(20*Mu2 + 4*mU32 + pow2(sbeta)*pow2(Yb)) +
        4*pow4(mQ3)) + (-mQ32 + mU32)*pow2(sbeta)*pow2(Yb)*pow4(Mu) + 2*(5*mQ32
        + mU32)*pow6(mD3)))/(mQ32*pow2(cbeta)*pow2(mD32 - Mu2)*pow3(mQ32 -
        mU32))) + (6*log(mU32/mQ32)*(-2*Xt*pow2(mD32 - Mu2)*((mD32 - mQ32 + Xb*
        Yt)*pow2(cbeta) + (mD32 - mQ32)*pow2(sbeta))*((-mA2 + mD32 + mU32)*Xb*
        pow2(cbeta) + (mD32 - mQ32)*Yb*pow2(sbeta))*pow3(mQ32 - mU32) - 2*(mD32
        - mQ32)*(2*mA2 - 2*mD32 - mQ32 - mU32)*(mQ32 - mU32)*Yb*pow2(mD32 -
        Mu2)*pow2(sbeta)*(Xb*Yt*pow2(cbeta) + (mD32 - mQ32)*(pow2(cbeta) +
        pow2(sbeta)))*pow3(Xt) + Xb2*pow2(cbeta)*pow2(mD32 - Mu2)*(2*mD32 +
        pow2(cbeta)*pow2(Yt))*pow4(mQ32 - mU32) - 2*cbeta*(mD32 - mQ32)*(mD32 -
        Mu2)*Xb*(cbeta*Mu2*Yt*(pow2(cbeta) + pow2(sbeta)) - mD32*(2*sbeta*std::abs(
        Mu) + cbeta*Yt*pow2(sbeta) + Yt*pow3(cbeta)))*pow4(mQ32 - mU32) + pow2(
        mD32 - mQ32)*(2*pow2(cbeta)*pow2(mD32 - Mu2)*pow2(sbeta) + pow2(mD32 -
        Mu2)*pow4(cbeta) + pow2(sbeta)*(-2*mD32*Mu2*(-1 + pow2(sbeta)) + (-1 +
        pow2(sbeta))*pow4(mD3) + pow2(sbeta)*pow4(Mu)))*pow4(mQ32 - mU32) +
        sbeta*Xt4*(4*cbeta*mD32*(mD32 - mQ32)*(mD32 - Mu2)*(mQ32 - mU32)*(mQ32
        + mU32)*Xb*std::abs(Mu) - 2*mD32*(mD32 - mU32)*(3*mD32 - 2*mQ32 - mU32)*
        sbeta*Xb2*pow2(cbeta)*pow2(mD32 - Mu2) - sbeta*pow2(mD32 - mQ32)*((-3*
        mA2 + mQ32 + 2*mU32)*pow2(sbeta)*pow2(Yb)*pow4(Mu) + pow4(mD3)*(-3*mA2*
        pow2(sbeta)*pow2(Yb) + 2*mU32*pow2(sbeta)*pow2(Yb) + mQ32*(-4*Mu2 +
        pow2(sbeta)*pow2(Yb)) - 2*Mu2*(4*mU32 + 3*pow2(sbeta)*pow2(Yb)) + pow4(
        mQ3) - pow4(mU3)) + mD32*Mu2*(2*mQ32*(Mu2 - pow2(sbeta)*pow2(Yb)) +
        Mu2*(4*mU32 + 3*pow2(sbeta)*pow2(Yb)) - 2*pow4(mQ3) + 2*(3*mA2*pow2(
        sbeta)*pow2(Yb) - 2*mU32*pow2(sbeta)*pow2(Yb) + pow4(mU3))) + (2*mQ32 +
        4*mU32 + 3*pow2(sbeta)*pow2(Yb))*pow6(mD3))) + Xt2*pow2(mQ32 - mU32)*(-
        (Xb2*pow2(cbeta)*pow2(mD32 - Mu2)*((-mA2 + mQ32)*pow2(cbeta)*pow2(Yt) +
        mD32*(2*mQ32 - 4*mU32*pow2(sbeta) + pow2(cbeta)*pow2(Yt)) + 4*pow2(
        sbeta)*pow4(mD3))) - 2*cbeta*(mD32 - mQ32)*(mD32 - Mu2)*Xb*(-(cbeta*
        Mu2*(-(mU32*Yb*pow2(sbeta)) + mQ32*Yt*(pow2(cbeta) + pow2(sbeta)) +
        mA2*(-(Yt*pow2(cbeta)) + (Yb - Yt)*pow2(sbeta)))) - mD32*(mU32*sbeta*(
        cbeta*sbeta*Yb + 4*std::abs(Mu)) + cbeta*(mA2 + Mu2)*(Yt*pow2(cbeta) + (-Yb
        + Yt)*pow2(sbeta)) - mQ32*(4*sbeta*std::abs(Mu) + cbeta*Yt*pow2(sbeta) + Yt*
        pow3(cbeta))) + (cbeta*(-Yb + Yt)*pow2(sbeta) + Yt*pow3(cbeta))*pow4(
        mD3)) + pow2(mD32 - mQ32)*(2*(mA2 - mD32 - mQ32)*pow2(cbeta)*pow2(mD32
        - Mu2)*pow2(sbeta) + (mA2 - mD32 - mQ32)*pow2(mD32 - Mu2)*pow4(cbeta) -
        pow2(sbeta)*(mD32*Mu2*(-4*mU32 - 2*mQ32*(-2 + pow2(sbeta)) + Mu2*(-2 +
        pow2(sbeta)) + 2*pow2(sbeta)*(mA2 + pow2(Yb))) - (-2*mU32 - mQ32*(-2 +
        pow2(sbeta)) + 2*Mu2*(-2 + pow2(sbeta)) + mA2*pow2(sbeta) + pow2(sbeta)
        *pow2(Yb))*pow4(mD3) - pow2(sbeta)*(mA2 - mQ32 + pow2(Yb))*pow4(Mu) + (
        -2 + pow2(sbeta))*pow6(mD3))))))/(pow2(cbeta)*pow2(mD32 - mQ32)*pow2(
        mD32 - Mu2)*pow2(sbeta)*pow4(mQ32 - mU32))) + (6*Xt4*((Xb2*(6*mQ32*
        dilog(1 - mD32/mQ32)*pow2(mD32 - mU32) - 6*mQ32*dilog(1 - mD32/mU32)*
        pow2(mD32 - mU32) + (mQ32 - mU32)*((-mD32 + mQ32)*(11*mQ32 + mU32) + 6*
        dilog(1 - mQ32/mU32)*(-(mQ32*mU32) + pow4(mQ3)))))/(-mD32 + mQ32) + (
        18*mD32*mQ32*Mu2*mU32 - 2*mQ32*(mD32 - Mu2)*(Mu2 - mU32)*(2*mQ32 - 3*
        Mu2 + mU32)*dilog(1 - Mu2/mQ32) + 8*mD32*mQ32*Mu2*mU32*dilog(1 - Mu2/
        mU32) + 6*mQ32*(mD32 - Mu2)*dilog(1 - mQ32/mU32)*pow2(cbeta)*pow2(mQ32
        - mU32) - 4*mA2*mD32*mQ32*mU32*pow2(sbeta) + 4*mA2*mQ32*Mu2*mU32*pow2(
        sbeta) - 10*mD32*mQ32*mU32*pow2(sbeta)*pow2(Yb) + 10*mQ32*Mu2*mU32*
        pow2(sbeta)*pow2(Yb) - 4*mQ32*mU32*pow4(mD3) - 21*mD32*Mu2*pow4(mQ3) +
        4*mD32*mU32*pow4(mQ3) + 4*mD32*Mu2*dilog(1 - Mu2/mU32)*pow4(mQ3) - 4*
        mD32*mU32*dilog(1 - Mu2/mU32)*pow4(mQ3) + 4*Mu2*mU32*dilog(1 - Mu2/
        mU32)*pow4(mQ3) - 48*mD32*mU32*pow2(cbeta)*pow4(mQ3) + 48*Mu2*mU32*
        pow2(cbeta)*pow4(mQ3) + 5*mA2*mD32*pow2(sbeta)*pow4(mQ3) - 5*mA2*Mu2*
        pow2(sbeta)*pow4(mQ3) - 2*mD32*mU32*pow2(sbeta)*pow4(mQ3) + 2*Mu2*mU32*
        pow2(sbeta)*pow4(mQ3) + 11*mD32*pow2(sbeta)*pow2(Yb)*pow4(mQ3) - 11*
        Mu2*pow2(sbeta)*pow2(Yb)*pow4(mQ3) + 5*pow4(mD3)*pow4(mQ3) - 14*mQ32*
        mU32*pow4(Mu) - 6*mD32*mQ32*dilog(1 - Mu2/mU32)*pow4(Mu) - 8*mQ32*mU32*
        dilog(1 - Mu2/mU32)*pow4(Mu) + 16*pow4(mQ3)*pow4(Mu) - 4*dilog(1 - Mu2/
        mU32)*pow4(mQ3)*pow4(Mu) - 8*mD32*mQ32*pow4(mU3) + 3*mD32*Mu2*pow4(mU3)
        + 6*mQ32*Mu2*pow4(mU3) - 2*mD32*mQ32*dilog(1 - Mu2/mU32)*pow4(mU3) + 2*
        mQ32*Mu2*dilog(1 - Mu2/mU32)*pow4(mU3) + 24*mD32*mQ32*pow2(cbeta)*pow4(
        mU3) - 24*mQ32*Mu2*pow2(cbeta)*pow4(mU3) - mA2*mD32*pow2(sbeta)*pow4(
        mU3) + mD32*mQ32*pow2(sbeta)*pow4(mU3) + mA2*Mu2*pow2(sbeta)*pow4(mU3)
        - mQ32*Mu2*pow2(sbeta)*pow4(mU3) - mD32*pow2(sbeta)*pow2(Yb)*pow4(mU3)
        + Mu2*pow2(sbeta)*pow2(Yb)*pow4(mU3) - pow4(mD3)*pow4(mU3) - 2*pow4(Mu)
        *pow4(mU3) + 4*mD32*pow6(mQ3) - 6*Mu2*pow6(mQ3) + 24*mD32*pow2(cbeta)*
        pow6(mQ3) - 24*Mu2*pow2(cbeta)*pow6(mQ3) + mD32*pow2(sbeta)*pow6(mQ3) -
        Mu2*pow2(sbeta)*pow6(mQ3) + 6*mQ32*dilog(1 - Mu2/mU32)*pow6(Mu))/((mD32
        - Mu2)*pow2(cbeta))))/(mQ32*pow4(mQ32 - mU32));
}

double ThresholdCalculator::getDeltaLambdaYtau4Yb2(int omitLogs) const
{
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double Xtau = p.Ae(2,2) - p.mu*p.vu/p.vd;
   const double Xtau2 = pow2(Xtau);
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mD3 = std::sqrt(mD32);
   const double mQ3 = std::sqrt(mQ32);
   const double mL3 = std::sqrt(p.ml2(2,2));
   const double mE3 = std::sqrt(p.me2(2,2));
   const double beta = std::atan(p.vu/p.vd);
   const double cbeta = std::cos(beta);
   const double sbeta = std::sin(beta);
   const double lmQ3MR = omitLogs*log(mQ32 / pow2(p.scale));
   const double eps = mQ3*0.01;

   if (std::abs(mE3 - mL3) < eps && std::abs(mD3 - mQ3) < eps) {
      return (-4*lmQ3MR*Xb*Xtau*(-6*pow2(mL3) + Xtau2))/(cbeta*sbeta*pow4(mL3));
   }

   if (std::abs(mE3 - mL3) < eps) {
      return (-4*Xb*Xtau*(-6*pow2(mL3) + Xtau2)*((-1 + lmQ3MR)*(mD32 - mQ32) + mD32*log(
        mD32/mQ32)))/(cbeta*(mD32 - mQ32)*sbeta*pow4(mL3));
   }

   if (std::abs(mD3 - mQ3) < eps) {
      const double mE32 = pow2(mE3);
      const double mL32 = pow2(mL3);
      return (-24*lmQ3MR*Xb*Xtau*(2*(-mE32 + mL32)*Xtau2 + log(mE32/mL32)*(mL32*Xtau2
        + mE32*(2*mL32 + Xtau2) - pow4(mE3) - pow4(mL3))))/(cbeta*sbeta*pow3(
        mE32 - mL32));
   }

   const double exact
      = (-24*Xb*Xtau*((-1 + lmQ3MR)*(mD32 - mQ32) + mD32*log(mD32/mQ32))*(2*
        Xtau2*(-pow2(mE3) + pow2(mL3)) + log(pow2(mE3)/pow2(mL3))*(Xtau2*pow2(
        mL3) + pow2(mE3)*(Xtau2 + 2*pow2(mL3)) - pow4(mE3) - pow4(mL3))))/(
        cbeta*(mD32 - mQ32)*sbeta*pow3(pow2(mE3) - pow2(mL3)));

   return exact;
}

double ThresholdCalculator::getDeltaLambdaYtau2Yb4(int omitLogs) const
{
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double Xtau = p.Ae(2,2) - p.mu*p.vu/p.vd;
   const double Xb2 = pow2(Xb);
   const double mL3 = std::sqrt(p.ml2(2,2));
   const double mE3 = std::sqrt(p.me2(2,2));
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mD3 = std::sqrt(mD32);
   const double mQ3 = std::sqrt(mQ32);
   const double lmE3MR = omitLogs*log(pow2(mE3)/pow2(p.scale));
   const double lmL3MR = omitLogs*log(pow2(mL3)/pow2(p.scale));
   const double lmD3MR = omitLogs*log(pow2(mD3)/pow2(p.scale));
   const double lmQ3MR = omitLogs*log(mQ32 / pow2(p.scale));
   const double beta = std::atan(p.vu/p.vd);
   const double cbeta = std::cos(beta);
   const double sbeta = std::sin(beta);
   const double eps = mQ3*0.01;

   if (std::abs(mE3 - mL3) < eps && std::abs(mD3 - mQ3) < eps) {
      return (-4*lmL3MR*Xb*(-6*mQ32 + Xb2)*Xtau)/(cbeta*sbeta*pow4(mQ3));
   }

   if (std::abs(mE3 - mL3) < eps) {
      return (24*lmL3MR*Xb*Xtau*(2*(mD32 - mQ32)*Xb2 + lmQ3MR*(mQ32*Xb2 + mD32*(2*mQ32
        + Xb2) - pow4(mD3) - pow4(mQ3)) + lmD3MR*(-(mQ32*Xb2) - mD32*(2*mQ32 +
        Xb2) + pow4(mD3) + pow4(mQ3))))/(cbeta*sbeta*pow3(mD32 - mQ32));
   }

   if (std::abs(mD3 - mQ3) < eps) {
      const double mE32 = pow2(mE3);
      const double mL32 = pow2(mL3);
      return (-4*(-mE32 + lmE3MR*mE32 + mL32 - lmL3MR*mL32)*Xb*(-6*mQ32 + Xb2)*Xtau)/(
        cbeta*(mE32 - mL32)*sbeta*pow4(mQ3));
   }

   const double exact
      = (24*Xb*Xtau*(-pow2(mE3) + lmE3MR*pow2(mE3) + pow2(mL3) -
        lmL3MR*pow2(mL3))*(2*(mD32 - mQ32)*Xb2 + lmQ3MR*(mQ32*
        Xb2 + mD32*(2*mQ32 + Xb2) - pow4(mD3) - pow4(mQ3)) + lmD3MR
        *(-(mQ32*Xb2) - mD32*(2*mQ32 + Xb2) + pow4(mD3) + pow4(mQ3))))/(cbeta*
        sbeta*(pow2(mE3) - pow2(mL3))*pow3(mD32 - mQ32));

   return exact;
}

//        at*as^n threshold corrections with limits and DR' -> MS shifts                //
/**
 * Returns delta g3_as in the MSbar scheme for a given mass limit
 * @param omitLogs an intiger to omit all log mu terms
 * @return delta g3_as in the MSbar scheme for a given mass limit
 */
double ThresholdCalculator::getDeltaG3Alphas(int omitLogs) const
{

   using std::log;

   const double Mst12 = pow2(p.MSt(0));
   const double m32 = pow2(p.MG);
   const double MR2 = pow2(p.scale);
   const double mU32 = p.mu2(2,2);
   const double mQ32 = p.mq2(2,2);
   const double lmsqMR = omitLogs * log(Mst12 / MR2) + log(msq2 / Mst12);
   const double lm3MR = omitLogs * log(Mst12 / MR2) + log(m32 / Mst12);
   const double lmU3MR = omitLogs * log(Mst12 / MR2) + log(mU32 / Mst12);
   const double lmQ3MR = omitLogs * log(Mst12 / MR2) + log(mQ32 / Mst12);

   return 1 / 2. - lm3MR - (lmQ3MR + 10 * lmsqMR + lmU3MR) / 12.;
}

/**
 * Returns delta yt_as in the MSbar scheme for a given mass limit
 * @param limit an integer key for a mass limit
 * @param omitLogs an integer key to omit all mu terms
 * @return delta yt_as in the MSbar scheme for a given mass limit
 */
double ThresholdCalculator::getDeltaYtAlphas(Limits limit, int omitLogs) const
{

   using std::log;

   const double Mst12 = pow2(p.MSt(0));
   const double m32 = pow2(p.MG);
   const double MR2 = pow2(p.scale);
   const double mU32 = p.mu2(2,2);
   const double mQ32 = p.mq2(2,2);
   const double m3 = std::sqrt(m32);
   const double mU3 = std::sqrt(mU32);
   const double mQ3 = std::sqrt(mQ32);
   const double Xt = p.Au(2,2) - p.mu * p.vd / p.vu;
   const double lmQ3MR = omitLogs * log(Mst12 / MR2) + log(mQ32 / Mst12);

   if (limit == Limits::DEGENERATE) limit = Limits::MQ3_EQ_MU3_EQ_M3;

   switch (limit) {
      case Limits::GENERAL:{
         return (2*(-2*lmQ3MR + (mQ32*mU32)/((-m32 + mQ32)*(m32 - mU32)) + pow4(m3)/((m32
        - mQ32)*(m32 - mU32)) - (mU32*log(mU32/pow2(mQ3))*(2*m32*(mQ32 - mU32)
        - mQ32*mU32 - 4*m3*mU32*Xt + 4*Xt*pow3(m3) + pow4(mU3)))/((-mQ32 +
        mU32)*pow2(m32 - mU32)) - (log(m32/pow2(mQ3))*pow3(m3)*(-4*mQ32*mU32*Xt
        + 4*m32*(mQ32 + mU32)*Xt - 2*(mQ32 + mU32)*pow3(m3) - 4*Xt*pow4(m3) +
        m3*(pow4(mQ3) + pow4(mU3)) + 2*pow5(m3)))/(pow2(m32 - mQ32)*pow2(m32 -
        mU32))))/3.;
      }
      case Limits::MQ3_EQ_MU3:{
         return (-2*((m32 - mQ32)*((-1 + 2*lmQ3MR)*m32 - (1 + 2*lmQ3MR)*mQ32 + 4*m3*Xt) +
        2*(m3 - 2*Xt)*log(m32/pow2(mQ3))*pow3(m3)))/(3.*pow2(m32 - mQ32));
      }
      case Limits::MQ3_EQ_M3:{
         return (-4*lmQ3MR - ((mQ32 - mU32)*(mQ32 - 3*mU32 - 8*mQ3*Xt) + 2*mU32*(-2*mQ32
        + mU32 - 4*mQ3*Xt)*log(mU32/pow2(mQ3)))/pow2(mQ32 - mU32))/3.;
      }
      case Limits::MU3_EQ_M3:{
         return (2*(-2*lmQ3MR + (-3*mQ32 + mU3*(mU3 - 8*Xt))/(2.*(mQ32 - mU32)) - (log(
        mU32/pow2(mQ3))*(-2*mQ32*mU3*(mU3 - 2*Xt) + pow4(mQ3) + 2*pow4(mU3)))/
        pow2(mQ32 - mU32)))/3.;
      }
      case Limits::MQ3_EQ_MU3_EQ_M3:{
         return (-4*(mQ3 + lmQ3MR*mQ3 - Xt))/(3.*mQ3);
      }
      default:
         break;
   };

   throw std::runtime_error("Mass limit not included!");
}

/**
  * Returns delta yt_as^2 in the MSbar scheme for a given mass limit
  * @param limit an integer key for a mass limit
  * @param omitLogs an integer key to omit all mu terms
  * @return delta yt_as^2 in the MSbar scheme for a given mass limit
  */
double ThresholdCalculator::getDeltaYtAlphas2(Limits limit, int omitLogs) const
{
   using std::log;
   using himalaya::dilog;

   const double Mst12 = pow2(p.MSt(0));
   const double m32 = pow2(p.MG);
   const double MR2 = pow2(p.scale);
   const double mU32 = p.mu2(2,2);
   const double mQ32 = p.mq2(2,2);
   const double m3 = std::sqrt(m32);
   const double mU3 = std::sqrt(mU32);
   const double mQ3 = std::sqrt(mQ32);
   const double msq = std::sqrt(msq2);
   const double Xt = p.Au(2,2) - p.mu * p.vd / p.vu;
   const double lmQ3MR = omitLogs * std::log(Mst12 / MR2) + std::log(mQ32 / Mst12);

   switch (limit) {
      case Limits::GENERAL: {
         const double logmqQ3 = std::log(msq2/mQ32);
         const double logU3Q3 = std::log(mU32/mQ32);
         const double logm3Q3 = std::log(m32/mQ32);

         return (4*pow2(lmQ3MR) - 60*pow2(logmqQ3) - (16*Xt*(-3*m3*(mQ32 + 10*msq2
        + mU32) + 14*pow3(m3)))/((m32 - mQ32)*(m32 - mU32)) + (60*(mQ32 - msq2)
        *dilog(1 - msq2/mQ32)*(-3*m32*(mQ32 - msq2) + mQ32*(mQ32 + msq2) + (8*
        m3*(m32 - mQ32)*(mQ32 - msq2)*Xt)/(mQ32 - mU32) - 2*pow4(m3)))/pow3(m32
        - mQ32) + (60*(msq2 - mU32)*dilog(1 - msq2/mU32)*(-3*m32*(msq2 - mU32)
        - mU32*(msq2 + mU32) + (8*m3*(m32 - mU32)*(msq2 - mU32)*Xt)/(-mQ32 +
        mU32) + 2*pow4(m3)))/pow3(m32 - mU32) + 60*(m32 - msq2)*dilog(1 - msq2/
        m32)*(1/(-m32 + mQ32) + 1/(-m32 + mU32) + (3*mU32)/pow2(m32 - mU32) - (
        8*m3*(mQ32*(msq2 - 2*mU32) + msq2*mU32 + m32*(mQ32 - 2*msq2 + mU32))*
        Xt)/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + (mQ32*(-3*m32 + 3*mQ32 + 4*
        msq2))/pow3(-m32 + mQ32) - msq2*(3*(1/pow2(m32 - mQ32) + 1/pow2(m32 -
        mU32)) + (4*mU32)/pow3(m32 - mU32)) + (4*pow4(mQ3))/pow3(m32 - mQ32) +
        (4*pow4(mU3))/pow3(m32 - mU32)) - (20*logmqQ3*(-24*m3*mQ32*msq2*
        mU32*Xt + 24*msq2*(mQ32 + mU32)*Xt*pow3(m3) + 3*mQ32*msq2*pow4(mU3) +
        pow4(mQ3)*(3*msq2*mU32 + 10*pow4(mU3)) + pow4(m3)*(-15*msq2*mU32 +
        mQ32*(-15*msq2 + 44*mU32) + 11*pow4(mQ3) + 11*pow4(mU3)) + 3*m32*((3*
        msq2 - 7*mU32)*pow4(mQ3) + 3*msq2*pow4(mU3) - mQ32*(4*msq2*mU32 + 7*
        pow4(mU3))) - 24*msq2*Xt*pow5(m3) + (-23*mQ32 + 18*msq2 - 23*mU32)*
        pow6(m3) + 12*pow8(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) + (-(mQ32*
        mU32*(6*mU32*(10*msq2 + mU32) + mQ32*(60*msq2 + 169*mU32) + 6*pow4(mQ3)
        )) + pow4(m3)*(4*mQ32*(75*msq2 - 28*mU32) + 300*msq2*mU32 + 93*pow4(
        mQ3) + 93*pow4(mU3)) - 6*(47*mQ32 + 60*msq2 + 47*mU32)*pow6(m3) - 2*
        m32*(2*(45*msq2 - 52*mU32)*pow4(mQ3) + 9*(10*msq2 + mU32)*pow4(mU3) -
        8*mQ32*(15*msq2*mU32 + 13*pow4(mU3)) + 9*pow6(mQ3)) + 291*pow8(m3))/(
        pow2(m32 - mQ32)*pow2(m32 - mU32)) - (2*(mQ32 - mU32)*dilog(1 - mU32/
        mQ32)*(-8*m3*mQ32*mU32*Xt*(-(mQ32*mU32) + 3*pow4(mQ3) + 3*pow4(mU3)) -
        16*Xt*(7*mQ32*mU32 + 4*pow4(mQ3) + 4*pow4(mU3))*pow5(m3) + (-76*mQ32*
        mU32 + 34*pow4(mQ3) + 34*pow4(mU3))*pow6(m3) + 7*pow4(mU3)*pow6(mQ3) +
        pow4(m3)*(65*mU32*pow4(mQ3) + 65*mQ32*pow4(mU3) - 29*pow6(mQ3) - 29*
        pow6(mU3)) + 7*pow4(mQ3)*pow6(mU3) + 8*Xt*pow3(m3)*(7*mU32*pow4(mQ3) +
        7*mQ32*pow4(mU3) + 3*pow6(mQ3) + 3*pow6(mU3)) + 80*(mQ32 + mU32)*Xt*
        pow7(m3) - 14*(mQ32 + mU32)*pow8(m3) + 3*mU32*pow8(mQ3) + 3*mQ32*pow8(
        mU3) + m32*(-34*pow4(mQ3)*pow4(mU3) - 26*mU32*pow6(mQ3) - 26*mQ32*pow6(
        mU3) + 9*pow8(mQ3) + 9*pow8(mU3)) - 40*Xt*pow9(m3) + 12*power10(m3)))/(
        pow3(m32 - mQ32)*pow3(m32 - mU32)) + (pow2(logU3Q3)*((128*m32*
        pow2(m32 - mU32)*pow2(Xt)*pow4(mU3))/pow2(mQ32 - mU32) + Xt*((8*m3*(
        mQ32 - mU32)*pow2(m32 - mU32)*(-(mQ32*mU32) - 5*m32*(mQ32 + mU32) + 5*
        pow4(m3) + 3*pow4(mQ3) + 3*pow4(mU3)))/pow2(m32 - mQ32) - (8*(m32 -
        mU32)*((69*mQ32*mU32 + pow4(mQ3) - 86*pow4(mU3))*pow5(m3) + m3*mU32*(
        10*mU32*pow4(mQ3) + 30*mU32*pow4(msq) - 60*msq2*pow4(mU3) + mQ32*(60*
        msq2*mU32 - 30*pow4(msq) + 4*pow4(mU3)) - 3*pow6(mQ3) - 27*pow6(mU3)) +
        pow3(m3)*(-11*mU32*pow4(mQ3) - 30*mU32*pow4(msq) + mQ32*(-60*msq2*mU32
        + 30*pow4(msq) - 59*pow4(mU3)) + 60*msq2*pow4(mU3) + 3*pow6(mQ3) + 99*
        pow6(mU3)) - 18*(mQ32 - mU32)*pow7(m3)))/pow2(mQ32 - mU32)) + ((mQ32 -
        mU32)*pow4(mU3)*(4*mQ32*mU32 + 3*pow4(mQ3) + 30*pow4(msq) + 21*pow4(
        mU3)) + 2*(mQ32*(30*msq2 - 34*mU32) - 30*msq2*mU32 + pow4(mQ3) - 351*
        pow4(mU3))*pow6(m3) - 2*m32*mU32*(18*mU32*pow4(mQ3) + 6*mU32*(-15*msq2*
        mU32 + 5*pow4(msq) + 4*pow4(mU3)) + mQ32*(90*msq2*mU32 - 30*pow4(msq) +
        25*pow4(mU3)) - 3*pow6(mQ3)) + pow4(m3)*(33*mU32*pow4(mQ3) + 90*mU32*
        pow4(msq) - 120*msq2*pow4(mU3) + mQ32*(120*msq2*mU32 - 90*pow4(msq) +
        141*pow4(mU3)) - 9*pow6(mQ3) + 347*pow6(mU3)) + (-44*mQ32 + 556*mU32)*
        pow8(m3) - 128*power10(m3))/(mQ32 - mU32) + ((m32 - mU32)*(-mQ32 +
        mU32)*((-76*mQ32*mU32 + 34*pow4(mQ3) + 34*pow4(mU3))*pow6(m3) + 7*pow4(
        mU3)*pow6(mQ3) + pow4(m3)*(65*mU32*pow4(mQ3) + 65*mQ32*pow4(mU3) - 29*
        pow6(mQ3) - 29*pow6(mU3)) + 7*pow4(mQ3)*pow6(mU3) - 14*(mQ32 + mU32)*
        pow8(m3) + 3*mU32*pow8(mQ3) + 3*mQ32*pow8(mU3) + m32*(-34*pow4(mQ3)*
        pow4(mU3) - 26*mU32*pow6(mQ3) - 26*mQ32*pow6(mU3) + 9*pow8(mQ3) + 9*
        pow8(mU3)) + 12*power10(m3)))/pow3(m32 - mQ32)))/pow4(m32 - mU32) + (2*
        logU3Q3*(10*logmqQ3*(-2*(3*msq2 + mU32)*pow4(m3) + 3*
        mU32*pow4(msq) + (4*m3*(m32 - mU32)*Xt*(m32*mU32 + 12*msq2*mU32 - 6*
        pow4(msq) - pow4(mU3)))/(-mQ32 + mU32) + 3*m32*(-6*msq2*mU32 + 3*pow4(
        msq) + pow4(mU3)) - pow6(mU3)) + (8*(m32 - mU32)*Xt*(-(mU32*(53*mQ32 +
        30*msq2 + 61*mU32)*pow3(m3)) + m3*mU32*(mQ32*(30*msq2 + 53*mU32) + 3*
        pow4(mQ3) + 3*pow4(mU3)) + (-16*mQ32 + 55*mU32)*pow5(m3) + 16*pow7(m3))
        )/((m32 - mQ32)*(mQ32 - mU32)) + (-(mQ32*(mQ32 - mU32)*(2*mQ32*(15*msq2
        - 64*mU32) + 3*pow4(mQ3) - 3*pow4(mU3))*pow4(mU3)) + mU32*pow4(m3)*(6*(
        30*msq2 + 157*mU32)*pow4(mQ3) - 5*mQ32*(42*msq2*mU32 - 5*pow4(mU3)) +
        30*msq2*pow4(mU3) + 163*pow6(mQ3) - 106*pow6(mU3)) + pow6(m3)*(-622*
        mU32*pow4(mQ3) - 5*mQ32*(18*msq2*mU32 + 197*pow4(mU3)) + 53*pow6(mQ3) +
        18*(5*msq2*pow4(mU3) + pow6(mU3))) + (764*mQ32*mU32 - 105*pow4(mQ3) +
        365*pow4(mU3))*pow8(m3) - m32*mU32*(-6*pow4(mQ3)*(25*msq2*mU32 - 41*
        pow4(mU3)) + (90*msq2 + 271*mU32)*pow6(mQ3) + 3*mQ32*(20*msq2*pow4(mU3)
        - 93*pow6(mU3)) + 9*pow8(mQ3) + 9*pow8(mU3)) + 4*(13*mQ32 - 77*mU32)*
        power10(m3))/((mQ32 - mU32)*pow2(m32 - mQ32))))/pow3(m32 - mU32) - 4*
        lmQ3MR*(20*logmqQ3 + (8*Xt*pow3(m3))/((m32 - mQ32)*(m32 - mU32))
        - (logU3Q3*(8*mU32*(mQ32 + 15*mU32)*pow4(m3) - 200*Xt*pow3(m3)*
        pow4(mU3) + (34*mQ32 - 98*mU32)*pow6(m3) + 23*(mQ32 - mU32)*pow6(mU3) +
        136*m3*Xt*pow6(mU3) + m32*(-69*mQ32*pow4(mU3) + 5*pow6(mU3)) + 64*Xt*
        pow7(m3)))/((-mQ32 + mU32)*pow3(m32 - mU32)) + (-463*m32*mQ32*mU32*(
        mQ32 + mU32) + 202*pow4(mQ3)*pow4(mU3) + pow4(m3)*(1044*mQ32*mU32 +
        257*pow4(mQ3) + 257*pow4(mU3)) - 573*(mQ32 + mU32)*pow6(m3) + 312*pow8(
        m3))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) - (logm3Q3*(-192*Xt*
        pow11(m3) + 154*pow12(m3) + 72*m32*(mQ32 + mU32)*pow4(mQ3)*pow4(mU3) -
        208*Xt*pow3(m3)*pow4(mQ3)*pow4(mU3) + 408*mQ32*mU32*(mQ32 + mU32)*Xt*
        pow5(m3) - 24*pow6(mQ3)*pow6(mU3) - pow6(m3)*(95*mU32*pow4(mQ3) + 95*
        mQ32*pow4(mU3) + 33*pow6(mQ3) + 33*pow6(mU3)) - pow4(m3)*(152*pow4(mQ3)
        *pow4(mU3) + 11*mU32*pow6(mQ3) + 11*mQ32*pow6(mU3)) - 200*Xt*(4*mQ32*
        mU32 + pow4(mQ3) + pow4(mU3))*pow7(m3) + (406*mQ32*mU32 + 163*pow4(mQ3)
        + 163*pow4(mU3))*pow8(m3) + 392*(mQ32 + mU32)*Xt*pow9(m3) - 288*(mQ32 +
        mU32)*power10(m3)))/(pow3(m32 - mQ32)*pow3(m32 - mU32))) + 2*pow2(logm3Q3)*
        (-15*(m32 - msq2)*(1/(m32 - mQ32) + 1/(m32 - mU32) - (3*
        mU32)/pow2(m32 - mU32) + (mQ32*(-3*m32 + 3*mQ32 + 4*msq2))/pow3(m32 -
        mQ32) + msq2*(3*(1/pow2(m32 - mQ32) + 1/pow2(m32 - mU32)) + (4*mU32)/
        pow3(m32 - mU32)) + (4*pow4(mQ3))/pow3(-m32 + mQ32) + (4*pow4(mU3))/
        pow3(-m32 + mU32)) + (64*pow2(Xt)*pow6(m3))/(pow2(m32 - mQ32)*pow2(m32
        - mU32)) + (8*m3*Xt*(-15*(m32 - mQ32)*(m32 - msq2)*(m32 - mU32)*(mQ32*(
        msq2 - 2*mU32) + msq2*mU32 + m32*(mQ32 - 2*msq2 + mU32)) - m32*(-142*
        m32*mQ32*mU32*(mQ32 + mU32) + 87*pow4(mQ3)*pow4(mU3) + pow4(m3)*(220*
        mQ32*mU32 + 57*pow4(mQ3) + 57*pow4(mU3)) - 82*(mQ32 + mU32)*pow6(m3) +
        27*pow8(m3))))/(pow3(m32 - mQ32)*pow3(m32 - mU32)) + (-878*(mQ32 +
        mU32)*pow14(m3) + 262*pow16(m3) - pow4(m3)*pow4(mQ3)*pow4(mU3)*(452*
        mQ32*mU32 + pow4(mQ3) + pow4(mU3)) + pow12(m3)*(2972*mQ32*mU32 + 885*
        pow4(mQ3) + 885*pow4(mU3)) + 144*m32*(mQ32 + mU32)*pow6(mQ3)*pow6(mU3)
        - 36*pow8(mQ3)*pow8(mU3) + pow8(m3)*(2404*pow4(mQ3)*pow4(mU3) + 1004*
        mU32*pow6(mQ3) + 1004*mQ32*pow6(mU3) + 49*pow8(mQ3) + 49*pow8(mU3)) -
        2*pow6(m3)*(184*pow4(mU3)*pow6(mQ3) + 184*pow4(mQ3)*pow6(mU3) + 79*
        mU32*pow8(mQ3) + 79*mQ32*pow8(mU3)) - 4*(733*mU32*pow4(mQ3) + 733*mQ32*
        pow4(mU3) + 80*pow6(mQ3) + 80*pow6(mU3))*power10(m3))/(pow4(m32 - mQ32)
        *pow4(m32 - mU32))) + (2*dilog(1 - m32/mQ32)*(-8*m3*(m32 - mU32)*Xt*(
        m32*(mQ32 - 37*mU32) - 7*mQ32*mU32 + 18*pow4(m3) + 3*pow4(mQ3) + 22*
        pow4(mU3)) + (128*pow12(m3) - mU32*pow4(mQ3)*(mU32*pow4(mQ3) + 19*mQ32*
        pow4(mU3) + 3*pow6(mQ3) - 23*pow6(mU3)) + pow6(m3)*(-251*mU32*pow4(mQ3)
        - 1025*mQ32*pow4(mU3) + 9*pow6(mQ3) - 13*pow6(mU3)) + (902*mQ32*mU32 +
        98*pow4(mQ3) + 280*pow4(mU3))*pow8(m3) + pow4(m3)*(417*pow4(mQ3)*pow4(
        mU3) - 151*mU32*pow6(mQ3) + 391*mQ32*pow6(mU3) + 20*pow8(mQ3) - 37*
        pow8(mU3)) - 2*(147*mQ32 + 173*mU32)*power10(m3) + m32*(41*pow4(mU3)*
        pow6(mQ3) - 167*pow4(mQ3)*pow6(mU3) + 41*mU32*pow8(mQ3) - 34*mQ32*pow8(
        mU3) - 9*power10(mQ3)))/pow2(m32 - mQ32)))/((mQ32 - mU32)*pow3(m32 -
        mU32)) + (2*dilog(1 - m32/mU32)*(8*m3*(m32 - mQ32)*Xt*(-7*mQ32*mU32 +
        m32*(-37*mQ32 + mU32) + 18*pow4(m3) + 22*pow4(mQ3) + 3*pow4(mU3)) + (-
        128*pow12(m3) + pow6(m3)*(1025*mU32*pow4(mQ3) + 251*mQ32*pow4(mU3) +
        13*pow6(mQ3) - 9*pow6(mU3)) + 19*pow6(mQ3)*pow6(mU3) - 2*(451*mQ32*mU32
        + 140*pow4(mQ3) + 49*pow4(mU3))*pow8(m3) - 23*pow4(mU3)*pow8(mQ3) +
        pow4(m3)*(-417*pow4(mQ3)*pow4(mU3) - 391*mU32*pow6(mQ3) + 151*mQ32*
        pow6(mU3) + 37*pow8(mQ3) - 20*pow8(mU3)) + pow4(mQ3)*pow8(mU3) + m32*
        mU32*(-41*pow4(mQ3)*pow4(mU3) + 167*mU32*pow6(mQ3) - 41*mQ32*pow6(mU3)
        + 34*pow8(mQ3) + 9*pow8(mU3)) + (346*mQ32 + 294*mU32)*power10(m3) + 3*
        mQ32*power10(mU3))/pow2(m32 - mU32)))/((mQ32 - mU32)*pow3(m32 - mQ32))
        + (2*logm3Q3*(8*(m32 - mQ32)*Xt*pow2(m32 - mU32)*pow3(m3)*(3*
        mU32*(10*msq2 + mU32) - 15*m32*(5*mQ32 + 4*msq2 + 5*mU32) + mQ32*(30*
        msq2 + 59*mU32) + 85*pow4(m3) + 3*pow4(mQ3)) - (m32 - mU32)*(864*pow12(
        m3) + 48*pow6(mQ3)*pow6(mU3) + 3*m32*mQ32*mU32*(5*(2*msq2 - 23*mU32)*
        pow4(mQ3) - 115*mQ32*pow4(mU3) + 10*msq2*pow4(mU3) + pow6(mQ3) + pow6(
        mU3)) - 9*pow6(m3)*((30*msq2 + 317*mU32)*pow4(mQ3) + 30*msq2*pow4(mU3)
        + mQ32*(-20*msq2*mU32 + 317*pow4(mU3)) + 39*pow6(mQ3) + 39*pow6(mU3)) +
        4*(60*msq2*mU32 + 2*mQ32*(30*msq2 + 511*mU32) + 365*pow4(mQ3) + 365*
        pow4(mU3))*pow8(m3) + pow4(m3)*(pow4(mQ3)*(-90*msq2*mU32 + 1894*pow4(
        mU3)) + (90*msq2 + 572*mU32)*pow6(mQ3) + 9*(10*msq2 + mU32)*pow6(mU3) +
        mQ32*(-90*msq2*pow4(mU3) + 572*pow6(mU3)) + 9*pow8(mQ3)) - 2*(971*mQ32
        + 90*msq2 + 971*mU32)*power10(m3)) + 10*(m32 - mU32)*logmqQ3*(-
        4*pow12(m3) + 42*mU32*pow4(mQ3)*pow6(m3) + 42*mQ32*pow4(mU3)*pow6(m3) -
        4*m3*(m32 - mQ32)*(m32 - mU32)*Xt*(-((mQ32 - 12*msq2 + mU32)*pow4(m3))
        + m32*(mQ32*mU32 - 12*pow4(msq)) + 6*mU32*pow4(msq) + 6*mQ32*(-2*msq2*
        mU32 + pow4(msq)) + pow6(m3)) - 14*mU32*pow4(m3)*pow6(mQ3) + 2*pow6(m3)
        *pow6(mQ3) - 14*mQ32*pow4(m3)*pow6(mU3) + 2*pow6(m3)*pow6(mU3) - 84*
        mQ32*mU32*pow8(m3) - 6*pow4(mQ3)*pow8(m3) - 6*pow4(mU3)*pow8(m3) + 20*
        mQ32*power10(m3) + 20*mU32*power10(m3) + 3*(m32 - msq2)*(-(mQ32*msq2*
        mU32*(pow4(mQ3) + pow4(mU3))) + (-8*msq2*mU32 + mQ32*(-8*msq2 + 30*
        mU32) + 3*pow4(mQ3) + 3*pow4(mU3))*pow6(m3) - pow4(m3)*((-9*msq2 + 15*
        mU32)*pow4(mQ3) - 9*msq2*pow4(mU3) + 3*mQ32*(2*msq2*mU32 + 5*pow4(mU3))
        + pow6(mQ3) + pow6(mU3)) + m32*(3*msq2*mU32*pow4(mQ3) + (-3*msq2 + 5*
        mU32)*pow6(mQ3) - 3*msq2*pow6(mU3) + mQ32*(3*msq2*pow4(mU3) + 5*pow6(
        mU3))) + (-8*mQ32 + 6*msq2 - 8*mU32)*pow8(m3) + 2*power10(m3))) + (2*
        logU3Q3*(604*mQ32*mU32*pow12(m3) - 204*mQ32*pow14(m3) - 180*
        mU32*pow14(m3) + 64*pow16(m3) + 64*mU32*pow2(m32 - mQ32)*pow2(m32 -
        mU32)*pow2(Xt)*pow4(m3) + 234*pow12(m3)*pow4(mQ3) + 122*pow12(m3)*pow4(
        mU3) - 239*pow4(mU3)*pow6(m3)*pow6(mQ3) - 13*pow4(mQ3)*pow6(m3)*pow6(
        mU3) - 53*pow4(m3)*pow6(mQ3)*pow6(mU3) + 503*pow4(mQ3)*pow4(mU3)*pow8(
        m3) + 401*mU32*pow6(mQ3)*pow8(m3) + 13*mQ32*pow6(mU3)*pow8(m3) + 43*
        pow4(m3)*pow4(mU3)*pow8(mQ3) - 76*mU32*pow6(m3)*pow8(mQ3) + 24*m32*
        pow6(mU3)*pow8(mQ3) + 13*pow8(m3)*pow8(mQ3) + 61*pow4(m3)*pow4(mQ3)*
        pow8(mU3) - 49*mQ32*pow6(m3)*pow8(mU3) - 6*m32*pow6(mQ3)*pow8(mU3) +
        30*pow8(m3)*pow8(mU3) - 6*pow8(mQ3)*pow8(mU3) - 2*(m32 - mQ32)*(m32 -
        mU32)*Xt*(68*pow11(m3) + 12*m3*pow4(mQ3)*pow6(mU3) + pow5(m3)*(-184*
        mU32*pow4(mQ3) - 215*mQ32*pow4(mU3) + 3*pow6(mU3)) + pow3(m3)*(105*
        pow4(mQ3)*pow4(mU3) - 17*mQ32*pow6(mU3)) + (353*mQ32*mU32 + 75*pow4(
        mQ3) + 120*pow4(mU3))*pow7(m3) - (137*mQ32 + 183*mU32)*pow9(m3)) - 779*
        mU32*pow4(mQ3)*power10(m3) - 369*mQ32*pow4(mU3)*power10(m3) - 101*pow6(
        mQ3)*power10(m3) - 31*pow6(mU3)*power10(m3) + 13*mQ32*pow4(m3)*power10(
        mU3) - 18*m32*pow4(mQ3)*power10(mU3) - 7*pow6(m3)*power10(mU3) + 6*
        pow6(mQ3)*power10(mU3)))/(mQ32 - mU32)))/(pow3(m32 - mQ32)*pow4(m32 -
        mU32)))/18.;
      }

      case Limits::MQ3_EQ_MU3: {
         const double logmqQ3 = std::log(msq2/mQ32);
         const double logm3Q3 = std::log(m32/mQ32);

         return (4*pow2(lmQ3MR) + (128*m3*Xt*(m32 - mQ32 + m3*Xt))/pow2(m32 - mQ32) - 60*
        pow2(logmqQ3) + (480*m3*msq2*(-mQ32 + msq2)*Xt*logmqQ3)/
        pow2(-(m32*mQ3) + pow3(mQ3)) + (32*(8*m3 - 9*Xt)*logm3Q3*pow3(m3)
        )/((m32 - mQ32)*mQ32) - (32*Xt*(-3*m3*(mQ32 + 5*msq2) + 7*pow3(m3)))/
        pow2(m32 - mQ32) + (60*(mQ32 - msq2)*dilog(1 - msq2/mQ32)*(-3*m32*(mQ32
        - msq2) + mQ32*(mQ32 + msq2) + (2*m3*(m32 - mQ32)*(mQ32 - msq2)*Xt)/
        mQ32 - 2*pow4(m3)))/pow3(m32 - mQ32) - (120*(m32 - msq2)*dilog(1 -
        msq2/m32)*(mQ32*msq2 + m32*(-5*mQ32 + 3*msq2) + 8*m3*(mQ32 - msq2)*Xt +
        pow4(m3)))/pow3(m32 - mQ32) + (18*m32*(mQ32 - 20*msq2) - 120*mQ32*msq2
        + 291*pow4(m3) - 181*pow4(mQ3))/pow2(m32 - mQ32) - (40*logmqQ3*(
        3*mQ32*msq2 + m32*(-11*mQ32 + 9*msq2) - 12*m3*msq2*Xt + 6*pow4(m3) + 5*
        pow4(mQ3)))/pow2(m32 - mQ32) + (60*(mQ32 - msq2)*dilog(1 - msq2/mQ32)*(
        2*m3*mQ32*(mQ32 - 9*msq2)*Xt + 2*(7*mQ32 + msq2)*Xt*pow3(m3) - 2*mQ32*
        pow4(m3) + (mQ32 + msq2)*pow4(mQ3) - 3*m32*(-(mQ32*msq2) + pow4(mQ3))))
        /(mQ32*pow3(m32 - mQ32)) - (16*m3*(55*m32*mQ32*Xt - 30*mQ32*msq2*Xt -
        32*mQ32*pow3(m3) + 16*Xt*pow4(m3) + 32*m3*pow4(mQ3) - 59*Xt*pow4(mQ3) -
        5*Xt*logmqQ3*(m32*mQ32 + 12*mQ32*msq2 - pow4(mQ3) - 6*pow4(msq))
        ))/pow2(-(m32*mQ3) + pow3(mQ3)) + (4*dilog(1 - m32/mQ32)*(-40*mQ32*Xt*
        pow3(m3) + 13*mQ32*pow4(m3) - 36*m32*pow4(mQ3) + 22*m3*Xt*pow4(mQ3) +
        18*Xt*pow5(m3) - 16*pow6(m3) + 15*pow6(mQ3)))/pow2(-(m32*mQ3) + pow3(
        mQ3)) + (4*dilog(1 - m32/mQ32)*(32*mQ32*Xt*pow3(m3) - 51*mQ32*pow4(m3)
        - 4*m32*pow4(mQ3) - 14*m3*Xt*pow4(mQ3) - 18*Xt*pow5(m3) + 16*pow6(m3) +
        15*pow6(mQ3)))/pow2(-(m32*mQ3) + pow3(mQ3)) - (8*lmQ3MR*((m32 - mQ32)*(
        -10*logmqQ3*pow2(-(m32*mQ3) + pow3(mQ3)) + 28*mQ32*Xt*pow3(m3) -
        188*mQ32*pow4(m3) + 293*m32*pow4(mQ3) - 68*m3*Xt*pow4(mQ3) + 32*Xt*
        pow5(m3) - 101*pow6(mQ3)) + mQ32*logm3Q3*(104*mQ32*Xt*pow3(m3) -
        57*mQ32*pow4(m3) - 36*m32*pow4(mQ3) - 96*Xt*pow5(m3) + 77*pow6(m3) +
        12*pow6(mQ3))))/(mQ32*pow3(-m32 + mQ32)) + (8*logm3Q3*(6*mQ32*(3*
        mQ32 + 20*msq2)*Xt*pow3(m3) + mQ32*(-291*mQ32 - 45*msq2 + 32*pow2(Xt))*
        pow4(m3) + 14*mQ32*Xt*pow5(m3) + 152*mQ32*pow6(m3) - 5*mQ32*logmqQ3*
        (24*m3*msq2*(-mQ32 + msq2)*Xt + 2*(mQ32 - 12*msq2)*Xt*pow3(m3) -
        (mQ32 - 6*msq2)*pow4(m3) + 9*m32*(2*mQ32*msq2 - pow4(msq)) - 3*mQ32*
        pow4(msq) - 2*Xt*pow5(m3) + pow6(m3)) - 12*m3*Xt*pow6(mQ3) + 15*m32*(-(
        msq2*pow4(mQ3)) + 9*pow6(mQ3)) - 68*Xt*pow7(m3) + 32*pow8(m3) - 12*
        pow8(mQ3)))/(mQ32*pow3(-m32 + mQ32)) + (4*pow2(logm3Q3)*(-120*m3*
        mQ32*(mQ32 - msq2)*msq2*Xt - 12*Xt*pow3(m3)*(19*pow4(mQ3) + 10*pow4(
        msq)) + pow4(m3)*(-60*mQ32*msq2 + 94*pow4(mQ3) + 45*pow4(msq)) + 40*(8*
        mQ32 + 3*msq2)*Xt*pow5(m3) + (-264*mQ32 - 30*msq2 + 32*pow2(Xt))*pow6(
        m3) + 6*m32*(15*msq2*pow4(mQ3) - 5*mQ32*pow4(msq) + 12*pow6(mQ3)) -
        108*Xt*pow7(m3) + 116*pow8(m3) - 3*(5*pow4(mQ3)*pow4(msq) + 6*pow8(mQ3)
        )))/pow4(m32 - mQ32))/18.;
      }

      case Limits::MQ3_EQ_M3: {
         const double logmqQ3 = std::log(msq2/mQ32);
         const double logU3Q3 = std::log(mU32/mQ32);
         const double logQ3U3 = -logU3Q3;

         return   (8064 + (384*(mQ32 + mU32 + 2*mQ3*Xt)*logQ3U3)/(mQ32 - mU32) +
        192*pow2(lmQ3MR) - 2880*pow2(logmqQ3) - (288*(mQ3 - 2*Xt)*(-2*
        mQ32 + 2*mU32 + (mQ32 + mU32)*logQ3U3))/(-(mQ3*mU32) + pow3(mQ3)
        ) - (2880*msq2*(-2*mQ32 + 2*msq2 + (-3*mQ32 + msq2)*logmqQ3)*(-(
        mQ3*mU32) - 2*mQ32*Xt + 2*msq2*Xt + pow3(mQ3)))/((mQ32 - msq2)*(mQ32 -
        mU32)*pow3(mQ3)) - (2880*(msq2 - mU32)*dilog(1 - msq2/mU32)*(3*mQ32*(
        msq2 - mU32) + mU32*(msq2 + mU32) + 8*mQ3*(msq2 - mU32)*Xt - 2*pow4(
        mQ3)))/pow3(mQ32 - mU32) + (180*(mQ32 - msq2)*(-3*mQ32 - msq2 - (8*mQ3*
        (mQ32 - msq2)*Xt)/(mQ32 - mU32))*dilog(1 - msq2/mQ32))/pow4(mQ3) - (
        192*Xt*(mQ32*(90*msq2 - 58*mU32) + 3*mU32*(10*msq2 + mU32) + 23*pow4(
        mQ3)))/(mQ3*pow2(mQ32 - mU32)) + (32*(mQ32*(-10*mQ32 + 30*msq2 + 3*mU32
        + 3*mQ32*logU3Q3) + 30*logmqQ3*(2*mQ32*msq2 - pow4(msq)))
        )/pow4(mQ3) + (96*(2*mQ32*mU32 - 3*pow4(mQ3) + logQ3U3*(3*mQ32*
        mU32 + pow4(mQ3)) + pow4(mU3)))/(-(mQ32*mU32) + pow4(mQ3)) - (12*(240*
        msq2*mU32 + 10*mQ32*(72*msq2 + 61*mU32) - 837*pow4(mQ3) + 19*pow4(mU3))
        )/pow2(mQ32 - mU32) - (240*logmqQ3*(72*mQ32*msq2*Xt + 24*msq2*
        mU32*Xt + 18*(2*msq2 - 5*mU32)*pow3(mQ3) + mQ3*(12*msq2*mU32 + 43*pow4(
        mU3)) + 47*pow5(mQ3)))/(mQ3*pow2(mQ32 - mU32)) + (24*(-(logU3Q3*
        pow3(mQ3)*(-34*mQ32*mU32 + 8*mQ3*mU32*Xt - 24*Xt*pow3(mQ3) + 13*pow4(
        mQ3) + 17*pow4(mU3))) + mQ3*(mQ32 - mU32)*(-(mQ32*(120*msq2 + 137*mU32)
        ) + 24*mQ3*(10*msq2 + mU32)*Xt - 104*Xt*pow3(mQ3) + 129*pow4(mQ3) + 12*
        (10*msq2*mU32 + pow4(mU3))) - 10*(mQ32 - mU32)*logmqQ3*(-12*mQ3*
        msq2*mU32 - 48*mQ32*msq2*Xt + (12*msq2 - mU32)*pow3(mQ3) + 24*Xt*pow4(
        msq) + pow5(mQ3))))/(pow2(mQ32 - mU32)*pow3(mQ3)) - (960*msq2*(-10*
        mQ32*msq2 + 7*pow4(mQ3) + 3*pow4(msq) + logmqQ3*(-3*mQ32*msq2 +
        6*pow4(mQ3) + pow4(msq))))/(-(msq2*pow4(mQ3)) + pow6(mQ3)) + (48*(10*(
        3*mQ32 - 3*mU32 - 8*mQ3*Xt)*logmqQ3*pow2(-(mQ3*mU32) + pow3(mQ3)
        ) - pow2(mQ32 - mU32)*(6*mU32*(10*msq2 + mU32) - 3*mQ32*(20*msq2 + 221*
        mU32) + 48*mQ3*(10*msq2 + mU32)*Xt - 1312*Xt*pow3(mQ3) + 1189*pow4(mQ3)
        ) + mQ32*logU3Q3*(144*mU32*Xt*pow3(mQ3) - 33*mU32*pow4(mQ3) -
        88*mQ3*Xt*pow4(mU3) + mQ32*(256*mU32*pow2(Xt) + 53*pow4(mU3)) + 8*Xt*
        pow5(mQ3) + 11*pow6(mQ3) - 15*pow6(mU3))))/(mQ32*pow3(mQ32 - mU32)) + (
        6*dilog(1 - mU32/mQ32)*(80*mU32*Xt*pow3(mQ3) + 31*mU32*pow4(mQ3) - 19*
        mQ32*pow4(mU3) - 24*mQ3*Xt*pow4(mU3) + 456*Xt*pow5(mQ3) + 119*pow6(mQ3)
        - 3*pow6(mU3)))/(-(mU32*pow4(mQ3)) + pow6(mQ3)) - (32*lmQ3MR*(120*logmqQ3*
        pow3(mQ32 - mU32) + (mQ32 - mU32)*(-2044*mQ32*mU32 - 1224*
        mQ3*mU32*Xt + 1176*Xt*pow3(mQ3) + 983*pow4(mQ3) + 1037*pow4(mU3)) + 6*
        logU3Q3*(64*mU32*Xt*pow3(mQ3) - 56*mU32*pow4(mQ3) - 5*mQ32*pow4(
        mU3) - 136*mQ3*Xt*pow4(mU3) + 64*Xt*pow5(mQ3) + 34*pow6(mQ3) + 23*pow6(
        mU3))))/pow3(mQ32 - mU32) + (24*logU3Q3*(-96*mQ32*mU32*(10*msq2
        + 19*mU32)*Xt + 1544*mU32*Xt*pow4(mQ3) - 3*pow3(mQ3)*(120*msq2*mU32 +
        59*pow4(mU3)) - 479*mU32*pow5(mQ3) + 512*Xt*pow6(mQ3) + 24*Xt*pow6(mU3)
        - 40*mQ3*logmqQ3*(4*mU32*Xt*pow3(mQ3) + 2*(3*msq2 + mU32)*pow4(
        mQ3) - 3*mU32*pow4(msq) - 3*mQ32*(-6*msq2*mU32 + 3*pow4(msq) + pow4(
        mU3)) - 4*mQ3*Xt*(-12*msq2*mU32 + 6*pow4(msq) + pow4(mU3)) + pow6(mU3))
        + mQ3*(-120*msq2*pow4(mU3) + 527*pow6(mU3)) + 209*pow7(mQ3)))/(mQ3*
        pow3(mQ32 - mU32)) + (4*(-2880*msq2*(msq2 + 2*mU32)*Xt*pow3(mQ3) +
        1440*mQ32*mU32*pow4(msq) + 2880*mQ3*mU32*Xt*pow4(msq) - 720*pow4(msq)*
        pow4(mU3) + pow4(mQ3)*(-720*pow4(msq) + 391*pow4(mU3)) + 480*(12*msq2 -
        5*mU32)*Xt*pow5(mQ3) + (-782*mU32 + 1536*pow2(Xt))*pow6(mQ3) + 2784*Xt*
        pow7(mQ3) + 487*pow8(mQ3)))/(pow2(mQ32 - mU32)*pow4(mQ3)) + (6*dilog(1
        - mQ32/mU32)*(-104*Xt*pow3(mQ3)*pow4(mU3) + 318*pow4(mQ3)*pow4(mU3) -
        2168*mU32*Xt*pow5(mQ3) + 1544*mU32*pow6(mQ3) + 16*mQ32*pow6(mU3) + 24*
        mQ3*Xt*pow6(mU3) + 2248*Xt*pow7(mQ3) - 2649*pow8(mQ3) + 3*pow8(mU3)))/(
        pow2(mQ32 - mU32)*pow4(mQ3)) - (180*dilog(1 - msq2/mQ32)*(16*msq2*mU32*
        (msq2 + mU32)*Xt*pow3(mQ3) - 8*mQ3*Xt*pow4(msq)*pow4(mU3) - 8*Xt*(-28*
        msq2*mU32 + 17*pow4(msq) + 17*pow4(mU3))*pow5(mQ3) + (90*msq2*mU32 -
        47*pow4(msq) + 39*pow4(mU3))*pow6(mQ3) + pow4(mQ3)*(-19*mU32*pow4(msq)
        + 6*msq2*pow4(mU3) - 13*pow6(mU3)) - pow4(msq)*pow6(mU3) + mQ32*(3*
        pow4(msq)*pow4(mU3) - 2*msq2*pow6(mU3)) + 16*(msq2 + mU32)*Xt*pow7(mQ3)
        + 17*(2*msq2 - 7*mU32)*pow8(mQ3) - 8*Xt*pow9(mQ3) + 29*power10(mQ3)))/(
        pow3(mQ32 - mU32)*pow4(mQ3)) - (3*pow2(logU3Q3)*(-2248*Xt*pow11(
        mQ3) + 2745*pow12(mQ3) - 3*pow12(mU3) + 240*mU32*Xt*(32*msq2*mU32 - 16*
        pow4(msq) + 15*pow4(mU3))*pow5(mQ3) - 4*pow6(mQ3)*(240*mU32*pow4(msq) -
        720*msq2*pow4(mU3) + 512*pow2(Xt)*pow4(mU3) + 215*pow6(mU3)) + 16*Xt*(-
        480*msq2*mU32 + 240*pow4(msq) - 677*pow4(mU3))*pow7(mQ3) + (-1920*msq2*
        mU32 + 1440*pow4(msq) + 3691*pow4(mU3))*pow8(mQ3) + 152*Xt*pow3(mQ3)*
        pow8(mU3) - pow4(mQ3)*(480*pow4(msq)*pow4(mU3) + 257*pow8(mU3)) + 8840*
        mU32*Xt*pow9(mQ3) - 2*(480*msq2 + 2621*mU32)*power10(mQ3) - 10*mQ32*
        power10(mU3) - 24*mQ3*Xt*power10(mU3)))/pow4(-(mQ3*mU32) + pow3(mQ3)))/
        864.;
      }

      case Limits::MU3_EQ_M3: {
         const double logmqQ3 = std::log(msq2/mQ32);
         const double logmqU3 = std::log(msq2/mU32);
         const double logU3Q3 = std::log(mU32/mQ32);

         return(32080 - (1536*(mQ32 + mU3*(mU3 + 2*Xt))*logU3Q3)/(mQ32 - mU32) +
        768*pow2(lmQ3MR) - 11520*pow2(logmqQ3) - (1152*(mU3 - 2*Xt)*(2*(
        mQ32 - mU32) + (mQ32 + mU32)*logU3Q3))/(-(mQ32*mU3) + pow3(mU3))
        - (11520*msq2*(2*(msq2 - mU32) + (msq2 - 3*mU32)*logmqU3)*(-(
        mQ32*mU3) + 2*msq2*Xt - 2*mU32*Xt + pow3(mU3)))/((-mQ32 + mU32)*(-msq2
        + mU32)*pow3(mU3)) - (48*(9*(80*msq2 - 93*mU32)*mU32 + 10*mQ32*(24*msq2
        + 61*mU32) + 19*pow4(mQ3)))/pow2(mQ32 - mU32) - (11520*(mQ32 - msq2)*
        dilog(1 - msq2/mQ32)*(msq2*mU3*(3*mU3 + 8*Xt) + mQ32*(msq2 - mU3*(3*mU3
        + 8*Xt)) + pow4(mQ3) - 2*pow4(mU3)))/pow3(mQ32 - mU32) + (720*(msq2 -
        mU32)*((mQ32 - mU32)*(msq2 + 3*mU32) + 8*mU3*(msq2 - mU32)*Xt)*dilog(1
        - msq2/mU32))/((mQ32 - mU32)*pow4(mU3)) - (768*Xt*(mQ32*(30*msq2 - 58*
        mU32) + 90*msq2*mU32 + 3*pow4(mQ3) + 23*pow4(mU3)))/(mU3*pow2(mQ32 -
        mU32)) + (384*(2*mQ32*mU32 + pow4(mQ3) - 3*pow4(mU3) + logU3Q3*(
        3*mQ32*mU32 + pow4(mU3))))/(-(mQ32*mU32) + pow4(mU3)) - (48*(2*mU3*(
        mQ32 - mU32)*(-120*msq2*mU3*(mU3 - 2*Xt) + mQ32*(120*msq2 + mU3*(-137*
        mU3 + 24*Xt)) + (129*mU3 - 104*Xt)*pow3(mU3) + 12*pow4(mQ3)) + logU3Q3*
        pow3(mU3)*(10*mQ32*mU3*(-47*mU3 + 32*Xt) + mU32*(227*mU32 -
        352*mU3*Xt - 256*pow2(Xt)) + 235*pow4(mQ3)) + 20*(mQ32 - mU32)*logmqQ3*
        (-12*msq2*mU32*(mU3 - 4*Xt) + mQ32*(12*msq2*mU3 + pow3(mU3))
        - 24*Xt*pow4(msq) - pow5(mU3))))/(pow2(mQ32 - mU32)*pow3(mU3)) - (960*
        logmqQ3*(36*msq2*mU32*(mU3 + 2*Xt) + 6*mQ32*(2*msq2*(mU3 + 2*Xt)
        - 15*pow3(mU3)) + 43*mU3*pow4(mQ3) + 47*pow5(mU3)))/(mU3*pow2(mQ32 -
        mU32)) - (128*(3 - logU3Q3)*((mQ32 - mU32)*(-60*msq2*mU32 + 30*
        pow4(msq) - 37*pow4(mU3)) + 8*Xt*pow5(mU3)))/((mQ32 - mU32)*pow4(mU3))
        + (48*(-2 + logU3Q3)*(3*(80*msq2*mU3 + 27*pow3(mU3))*pow4(mQ3) +
        480*mU32*Xt*pow4(msq) + 240*msq2*(mU3 - 4*Xt)*pow4(mU3) - 2*mQ32*(240*
        msq2*mU32*(mU3 - 2*Xt) + 240*Xt*pow4(msq) + (81*mU3 - 184*Xt)*pow4(mU3)
        ) + (65*mU32 - 432*mU3*Xt - 256*pow2(Xt))*pow5(mU3)))/(pow2(mQ32 -
        mU32)*pow3(mU3)) - (96*(20*mU32*(3*mQ32 + mU3*(-3*mU3 + 8*Xt))*logmqQ3*
        pow2(mQ32 - mU32) + 2*pow2(mQ32 - mU32)*(-60*msq2*mU3*(mU3 -
        8*Xt) + mQ32*(60*msq2 - 663*mU32 + 48*mU3*Xt) + 41*(29*mU3 - 32*Xt)*
        pow3(mU3) + 6*pow4(mQ3)) + mU32*logU3Q3*(mQ32*mU32*(-1177*mU32 +
        864*mU3*Xt - 768*pow2(Xt)) + mU3*(625*mU3 - 656*Xt)*pow4(mQ3) + (563*
        mU32 - 336*mU3*Xt + 256*pow2(Xt))*pow4(mU3) - 43*pow6(mQ3))))/(mU32*
        pow3(-mQ32 + mU32)) + (96*logU3Q3*(mQ32*mU32*(3971*mU32 - 3328*
        mU3*Xt + 1280*pow2(Xt)) + mU3*(-2867*mU3 + 1856*Xt)*pow4(mQ3) + (-1665*
        mU32 + 1728*mU3*Xt - 256*pow2(Xt))*pow4(mU3) + 625*pow6(mQ3)))/pow3(
        mQ32 - mU32) - (128*lmQ3MR*(120*logmqQ3*pow3(mQ32 - mU32) + (
        mQ32 - mU32)*(-4*mQ32*mU3*(511*mU3 + 306*Xt) + (983*mU3 + 1176*Xt)*
        pow3(mU3) + 1037*pow4(mQ3)) + 6*logU3Q3*(2*mQ32*(-61*mU3 + 32*
        Xt)*pow3(mU3) + mU3*(61*mU3 - 136*Xt)*pow4(mQ3) + 8*(7*mU3 + 8*Xt)*
        pow5(mU3) + pow6(mQ3))))/pow3(mQ32 - mU32) + (128*(-((3*mQ32 + 30*msq2
        - 10*mU32)*(mQ32 - mU32)*mU32) + 30*msq2*(msq2 - 2*mU32)*(mQ32 - mU32)*
        logmqQ3 + 8*(-4*mQ32 + mU3*(4*mU3 + Xt))*logU3Q3*pow4(
        mU3)))/(-(mQ32*pow4(mU3)) + pow6(mU3)) + (24*dilog(1 - mU32/mQ32)*(-(
        mQ32*(31*mU3 + 80*Xt)*pow3(mU3)) + mU3*(19*mU3 + 24*Xt)*pow4(mQ3) - (
        119*mU3 + 456*Xt)*pow5(mU3) + 3*pow6(mQ3)))/(-(mQ32*pow4(mU3)) + pow6(
        mU3)) - (3840*msq2*(-10*msq2*mU32 + 3*pow4(msq) + 7*pow4(mU3) + logmqU3*
        (-3*msq2*mU32 + pow4(msq) + 6*pow4(mU3))))/(-(msq2*pow4(mU3)
        ) + pow6(mU3)) + (3*pow2(logU3Q3)*(120*(mU3 - 8*Xt)*pow3(mU3)*
        pow4(msq) + pow4(mQ3)*(240*msq2*mU32 + 120*pow4(msq) - 1267*pow4(mU3))
        + 240*msq2*(mU3 + 8*Xt)*pow5(mU3) - 2*mQ32*(240*msq2*(mU3 + 4*Xt)*pow3(
        mU3) + 120*mU3*(mU3 - 4*Xt)*pow4(msq) + (2573*mU3 - 4832*Xt)*pow5(mU3))
        + (6413*mU32 - 7616*mU3*Xt - 512*pow2(Xt))*pow6(mU3)))/(pow2(mQ32 -
        mU32)*pow4(mU3)) + (24*dilog(1 - mU32/mQ32)*(2*(159*mU3 - 52*Xt)*pow3(
        mU3)*pow4(mQ3) + 8*mQ32*(193*mU3 - 271*Xt)*pow5(mU3) + 8*mU3*(2*mU3 +
        3*Xt)*pow6(mQ3) + (-2649*mU3 + 2248*Xt)*pow7(mU3) + 3*pow8(mQ3)))/(
        pow2(mQ32 - mU32)*pow4(mU3)) + (24*logU3Q3*(-10*logmqQ3*
        pow2(mQ32 - mU32)*((mQ32 - mU32)*(6*msq2*mU32 + 3*pow4(msq) - 5*pow4(
        mU3)) + 8*mU3*Xt*(-6*msq2*mU32 + 3*pow4(msq) + 2*pow4(mU3))) + 8*Xt*
        pow3(mU3)*(5*(6*msq2 + 113*mU32)*pow4(mQ3) + 30*msq2*pow4(mU3) - 3*
        mQ32*(20*msq2*mU32 + 409*pow4(mU3)) + 3*pow6(mQ3) + 627*pow6(mU3)) +
        mU32*(-10*pow4(mQ3)*(9*msq2*mU32 + 325*pow4(mU3)) + 6*(5*msq2 + 2*mU32)
        *pow6(mQ3) - 3*(10*msq2 + 971*mU32)*pow6(mU3) + 18*mQ32*(5*msq2*pow4(
        mU3) + 338*pow6(mU3)) + 3*pow8(mQ3))))/(pow3(mQ32 - mU32)*pow4(mU3)) +
        (3*pow2(logU3Q3)*(4*(240*msq2*(mU3 + 2*Xt)*pow3(mU3) + 120*mU3*(
        mU3 - 2*Xt)*pow4(msq) + (-8781*mU3 + 12800*Xt)*pow5(mU3))*pow6(mQ3) +
        pow4(mQ3)*(-240*(11*mU3 - 12*Xt)*pow3(mU3)*pow4(msq) + 480*msq2*(21*mU3
        + 52*Xt)*pow5(mU3) + 2*(49527*mU32 - 43840*mU3*Xt + 9984*pow2(Xt))*
        pow6(mU3)) + (-240*msq2*mU32*(17*mU3 + 8*Xt) + (23053*mU32 - 7808*mU3*
        Xt - 512*pow2(Xt))*pow3(mU3) + 120*(47*mU3 + 136*Xt)*pow4(msq))*pow7(
        mU3) + (-240*msq2*mU32 - 120*pow4(msq) + 1837*pow4(mU3))*pow8(mQ3) -
        20*mQ32*(24*(7*mU3 + 38*Xt)*pow4(msq)*pow5(mU3) + 48*msq2*(7*mU3 + 26*
        Xt)*pow7(mU3) + (4377*mU32 - 2624*mU3*Xt - 256*pow2(Xt))*pow8(mU3))))/
        pow4(-(mQ32*mU3) + pow3(mU3)) + (720*dilog(1 - msq2/mU32)*((47*mU3 +
        136*Xt)*pow4(msq)*pow5(mU3) - pow4(mQ3)*(2*msq2*(3*mU3 + 8*Xt)*pow3(
        mU3) + mU3*(3*mU3 - 8*Xt)*pow4(msq) + (39*mU3 - 136*Xt)*pow5(mU3)) + (
        2*msq2*mU32 + pow4(msq) + 13*pow4(mU3))*pow6(mQ3) - 2*msq2*(17*mU3 + 8*
        Xt)*pow7(mU3) + mQ32*((19*mU3 - 16*Xt)*pow3(mU3)*pow4(msq) - 2*msq2*(
        45*mU3 + 112*Xt)*pow5(mU3) + (119*mU3 - 16*Xt)*pow7(mU3)) + (-29*mU3 +
        8*Xt)*pow9(mU3)))/(pow3(-mQ32 + mU32)*pow4(mU3)) - (6*logU3Q3*(
        4*(mQ32 - mU32)*mU32*(-2*(45*msq2*mU3*(mU3 + 8*Xt) + (11243*mU3 - 5748*
        Xt)*pow3(mU3))*pow4(mQ3) + 30*msq2*(47*mU3 + 104*Xt)*pow5(mU3) + 2*
        mQ32*(15*msq2*(19*mU3 + 48*Xt)*pow3(mU3) + 4*(3693*mU3 - 1951*Xt)*pow5(
        mU3)) + (30*msq2 + 8*mU3*(497*mU3 - 9*Xt))*pow6(mQ3) + 3*(-3807*mU3 +
        968*Xt)*pow7(mU3) + 3*pow8(mQ3)) + logU3Q3*pow4(mU3)*(2*mU32*(
        25391*mU32 - 28144*mU3*Xt + 2816*pow2(Xt))*pow4(mQ3) + 4*mQ32*(-12005*
        mU32 + 14088*mU3*Xt + 768*pow2(Xt))*pow4(mU3) - 4*mU3*(4621*mU3 - 5624*
        Xt)*pow6(mQ3) + (15341*mU32 - 19488*mU3*Xt - 512*pow2(Xt))*pow6(mU3) +
        1149*pow8(mQ3)) - 40*(mQ32 - mU32)*logmqQ3*(3*(47*mU3 + 136*Xt)*
        pow4(msq)*pow5(mU3) - pow4(mQ3)*(6*msq2*(3*mU3 + 8*Xt)*pow3(mU3) + 3*
        mU3*(3*mU3 - 8*Xt)*pow4(msq) + (33*mU3 - 80*Xt)*pow5(mU3)) + (6*msq2*
        mU32 + 3*pow4(msq) + 11*pow4(mU3))*pow6(mQ3) - 6*msq2*(17*mU3 + 8*Xt)*
        pow7(mU3) + mQ32*((57*mU3 - 48*Xt)*pow3(mU3)*pow4(msq) - 6*msq2*(45*mU3
        + 112*Xt)*pow5(mU3) + (49*mU3 - 96*Xt)*pow7(mU3)) + (-27*mU3 + 16*Xt)*
        pow9(mU3))))/pow4(-(mQ32*mU3) + pow3(mU3)))/3456.;
      }

      case Limits::MQ3_EQ_MU3_EQ_M3: {
         const double logmqQ3 = std::log(msq2/mQ32);

         return (1835*mQ32*msq2 - 232*mQ3*msq2*Xt + 780*mQ32*msq2*logmqQ3 - 360*
        mQ3*msq2*Xt*logmqQ3 - 4*lmQ3MR*mQ3*(mQ32 - msq2)*(335*mQ3 + 104*
        Xt + 60*mQ3*logmqQ3) + 96*mQ32*pow2(Xt) - 96*msq2*pow2(Xt) +
        180*mQ32*msq2*pow2(logmqQ3) + 232*Xt*pow3(mQ3) - 120*Xt*logmqQ3*pow3(mQ3) -
        1835*pow4(mQ3) - 540*logmqQ3*pow4(mQ3) -
        180*pow2(logmqQ3)*pow4(mQ3) - 360*dilog(1 - msq2/mQ32)*(-(mQ32*
        msq2) + pow4(mQ3)) + 12*pow2(lmQ3MR)*(-(mQ32*msq2) + pow4(mQ3)))/(54.*
        mQ32*(mQ32 - msq2));
      }

      case Limits::DEGENERATE: {
         return (-2075*mQ32 + 712*mQ3*Xt - 4*lmQ3MR*mQ3*(335*mQ3 + 104*Xt) + 12*mQ32*
        pow2(lmQ3MR) + 96*pow2(Xt))/(54.*mQ32);
      }

      default:
         break;
   };

   throw std::runtime_error("Mass limit not included!");
}

/**
 * Returns delta lambda_at in the MSbar scheme for a given mass limit
 * @param limit an integer key for a mass limit
 * @param omitLogs an integer key to omit all mu terms
 * @return delta lambda_at in the MSbar scheme for a given mass limit
 */
double ThresholdCalculator::getDeltaLambdaAlphat(Limits limit, int omitLogs) const
{
   const double Mst12 = pow2(p.MSt(0));
   const double MR2 = pow2(p.scale);
   const double mQ32 = p.mq2(2,2);
   const double Xt = p.Au(2,2) - p.mu * p.vd / p.vu;
   const double lmQ3MR = omitLogs * std::log(Mst12 / MR2) + std::log(mQ32 / Mst12);

   // Translate limits to occurring ones
   if (limit == Limits::DEGENERATE)
      limit = Limits::MQ3_EQ_MU3_EQ_M3;
   if (limit == Limits::MQ3_EQ_M3 || limit == Limits::MU3_EQ_M3)
      limit = Limits::GENERAL;
   if (limit == Limits::MQ3_EQ_MU3_EQ_M3)
      limit = Limits::MQ3_EQ_MU3;

   switch (limit) {
      case Limits::GENERAL: {
         const double mU32 = p.mu2(2,2);
         return 12*pow4(Xt)/pow2(mQ32 - mU32) + 12*lmQ3MR + (6 - 12*pow2(Xt)/(mQ32 - mU32)
            + (6*(mQ32 + mU32)/pow3(mQ32 - mU32))*pow4(Xt))*std::log(mU32/mQ32);
      }
      case Limits::MQ3_EQ_MU3:
         return 12*lmQ3MR + 12*pow2(Xt)/mQ32 - pow4(Xt)/pow2(mQ32);
      default:
         break;
   };

   throw std::runtime_error("Mass limit not included!");
}

/**
 * Returns delta lambda_atas in the MSbar scheme for a given mass limit
 * @param limit an integer key for a mass limit
 * @param omitLogs an integer key to omit all mu terms
 * @return delta lambda_atas in the MSbar scheme for a given mass limit
 */
double ThresholdCalculator::getDeltaLambdaAlphatAlphas(Limits limit, int omitLogs) const
{
   using std::log;
   using himalaya::dilog;

   const double Mst12 = pow2(p.MSt(0));
   const double m32 = pow2(p.MG);
   const double MR2 = pow2(p.scale);
   const double mU32 = p.mu2(2,2);
   const double mQ32 = p.mq2(2,2);
   const double m3 = std::sqrt(m32);
   const double mU3 = std::sqrt(mU32);
   const double mQ3 = std::sqrt(mQ32);
   const double Xt = p.Au(2,2) - p.mu * p.vd / p.vu;
   const double lmQ3MR = omitLogs * log(Mst12 / MR2) + log(mQ32 / Mst12);

   if (limit == Limits::DEGENERATE)
      limit = Limits::MQ3_EQ_MU3_EQ_M3;

   switch (limit) {
      case Limits::GENERAL:{
         return 16*(-6*pow2(lmQ3MR) - (4*m32*pow2(Xt))/(pow2(mQ3)*pow2(mU3)) - (32*m3*
        pow3(Xt))/pow2(mQ32 - mU32) - (2*dilog(1 - m32/pow2(mQ3))*(-8*(m32 -
        mQ32)*Xt*pow2(mQ32 - mU32)*pow3(m3) - 4*m3*(2*m32 - mQ32 - mU32)*pow2(
        m32 - mQ32)*pow3(Xt) + 2*pow3(mQ32 - mU32)*pow4(m3) - (2*m32 - mQ32 -
        mU32)*pow2(m32 - mQ32)*pow4(Xt)))/(pow2(m32 - mQ32)*pow3(mQ32 - mU32))
        - (2*dilog(1 - m32/pow2(mU3))*(8*(m32 - mU32)*Xt*pow2(mQ32 - mU32)*
        pow3(m3) + 4*m3*(2*m32 - mQ32 - mU32)*pow2(m32 - mU32)*pow3(Xt) + 2*
        pow3(mQ32 - mU32)*pow4(m3) + (2*m32 - mQ32 - mU32)*pow2(m32 - mU32)*
        pow4(Xt)))/(pow2(m32 - mU32)*pow3(mQ32 - mU32)) + (-6*m32*mQ32*mU32*(
        mQ32 + mU32) + 3*pow4(mQ3)*pow4(mU3) + pow4(m3)*(9*mQ32*mU32 + 2*pow4(
        mQ3) + 2*pow4(mU3)) - 2*(mQ32 + mU32)*pow6(m3))/((-m32 + mQ32)*(m32 -
        mU32)*pow2(mQ3)*pow2(mU3)) + (2*pow4(Xt)*(-5*m32*mQ32*mU32*(mQ32 +
        mU32) + 5*pow4(mQ3)*pow4(mU3) - pow4(m3)*(-5*mQ32*mU32 + pow4(mQ3) +
        pow4(mU3)) + (mQ32 + mU32)*pow6(m3)))/((m32 - mQ32)*(m32 - mU32)*pow2(
        mQ3)*pow2(mU3)*pow2(mQ32 - mU32)) + (2*lmQ3MR*(((mQ32 - mU32)*(8*m3*
        mQ32*mU32*pow3(Xt) + 2*mQ32*mU32*(-2*mQ32*mU32 + pow4(mQ3) + pow4(mU3)
        - 3*pow4(Xt)) - m32*(pow2(-(mU3*pow2(Xt)) + pow3(mU3)) - (mU32 + 2*
        pow2(Xt))*pow4(mQ3) + mQ32*(4*mU32*pow2(Xt) - pow4(mU3) + pow4(Xt)) +
        pow6(mQ3))))/(pow2(mQ3)*pow2(mU3)) - log(mU32/pow2(mQ3))*(mU32*(-4*m3 +
        3*Xt)*pow3(Xt) + (-9*mU32 + 4*m3*Xt - 6*pow2(Xt))*pow4(mQ3) + 2*(2*m3 -
        3*Xt)*Xt*pow4(mU3) + mQ32*(4*mU32*Xt*(-2*m3 + 3*Xt) + (-4*m3 + 3*Xt)*
        pow3(Xt) + 9*pow4(mU3)) + 2*m32*pow4(Xt) + 3*pow6(mQ3) - 3*pow6(mU3))))
        /pow3(mQ32 - mU32) - (pow2(log(mU32/pow2(mQ3)))*(4*(m32 - mU32)*Xt*(-(
        m3*mU32) + 2*pow3(m3))*pow3(mQ32 - mU32) + 4*m3*(m32 - mU32)*(mQ32 -
        mU32)*pow3(Xt)*(3*mQ32*mU32 - m32*(mQ32 + 3*mU32) + 2*pow4(m3) - pow4(
        mU3)) - 2*pow2(mQ32 - mU32)*pow2(Xt)*((mQ32 - 3*mU32)*pow4(m3) + 2*(
        mQ32 - 2*mU32)*pow4(mU3) + m32*(-4*mQ32*mU32 + 8*pow4(mU3))) + (-4*m32*
        mU32 + 3*pow4(m3) + 2*pow4(mU3))*pow4(mQ32 - mU32) - 4*m3*(m32 - mU32)*
        mU32*(mQ32 + mU32)*pow5(Xt) + pow4(Xt)*(-8*mQ32*mU32*pow4(m3) + (-4*
        mQ32*mU32 + pow4(mQ3) - 5*pow4(mU3))*pow4(mU3) + 2*(mQ32 - mU32)*pow6(
        m3) + m32*(-2*mU32*pow4(mQ3) + 10*mQ32*pow4(mU3) + 8*pow6(mU3)))))/(
        pow2(m32 - mU32)*pow4(mQ32 - mU32)) + (log(mU32/pow2(mQ3))*(8*m3*(m32 -
        mQ32)*Xt*pow2(m32 - mU32)*pow2(mQ32 - mU32) - 8*m3*(m32 - mQ32)*(mQ32 +
        3*mU32)*pow2(m32 - mU32)*pow3(Xt) - 2*(m32 - mU32)*pow2(mQ32 - mU32)*
        pow2(Xt)*(3*mQ32*mU32 - 2*m32*(3*mQ32 + 2*mU32) + 7*pow4(m3)) + 8*m3*(
        m32 - mQ32)*(m32 - mU32)*mU32*pow5(Xt) + pow3(mQ32 - mU32)*(-((6*mQ32 +
        7*mU32)*pow4(m3)) - 2*mQ32*pow4(mU3) + m32*(5*mQ32*mU32 + 3*pow4(mU3))
        + 7*pow6(m3)) + pow4(Xt)*(-3*mQ32*(mQ32 + 5*mU32)*pow4(mU3) - pow4(m3)*
        (15*mQ32*mU32 + 6*pow4(mQ3) + 29*pow4(mU3)) + (3*mQ32 + 7*mU32)*pow6(
        m3) + m32*(7*mU32*pow4(mQ3) + 31*mQ32*pow4(mU3) + 16*pow6(mU3)) + 4*
        pow8(m3))))/((m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) - (m32*
        log(m32/pow2(mQ3))*(-4*m32*(m32 - mQ32)*(m32 - mU32)*(m32 - mQ32 -
        mU32)*pow2(Xt)*pow3(mQ32 - mU32) - 8*m3*(m32 - mQ32)*mQ32*(m32 - mU32)*
        (mQ32 - mU32)*mU32*pow5(Xt) + m3*mQ32*mU32*log(mU32/pow2(mQ3))*(-4*(m32
        - mQ32)*(m32 - mU32)*(2*m32 - mQ32 - mU32)*Xt*pow2(mQ32 - mU32) - 8*(
        m32 - mQ32)*(m32 - mU32)*pow3(Xt)*(4*mQ32*mU32 - 2*m32*(mQ32 + mU32) +
        2*pow4(m3) - pow4(mQ3) - pow4(mU3)) - 2*m3*pow2(mQ32 - mU32)*pow2(Xt)*(
        -2*m32*(mQ32 + mU32) + 2*pow4(m3) + pow4(mQ3) + pow4(mU3)) + m3*(2*m32
        - mQ32 - mU32)*pow4(mQ32 - mU32) + m3*(mQ32 + mU32)*(-2*m32*(mQ32 +
        mU32) + 2*pow4(m3) + pow4(mQ3) + pow4(mU3))*pow4(Xt) - 4*(m32 - mQ32)*(
        m32 - mU32)*(mQ32 + mU32)*pow5(Xt)) + m32*pow3(mQ32 - mU32)*(2*m32*(
        mQ32 + mU32)*pow2(mQ32 - mU32) + pow4(m3)*(2*mQ32*mU32 - 4*pow4(mQ3) -
        4*pow4(mU3)) + mQ32*mU32*(pow4(mQ3) + pow4(mU3)) + 2*(mQ32 + mU32)*
        pow6(m3)) + 2*(mQ32 - mU32)*pow4(Xt)*(pow3(mQ32 + mU32)*pow4(m3) + (
        mQ32 + mU32)*pow4(mQ3)*pow4(mU3) - m32*mQ32*mU32*(4*mQ32*mU32 + pow4(
        mQ3) + pow4(mU3)) - 2*(mQ32*mU32 + pow4(mQ3) + pow4(mU3))*pow6(m3) + (
        mQ32 + mU32)*pow8(m3))))/(pow2(mQ3)*pow2(m32 - mQ32)*pow2(mU3)*pow2(m32
        - mU32)*pow3(mQ32 - mU32)));
      }
      case Limits::MQ3_EQ_MU3:{
         return (-8*(48*(m3 - 2*Xt)*dilog(1 - m32/pow2(mQ3))*pow3(m3)*pow6(mQ3) + 4*log(
        m32/pow2(mQ3))*pow3(m3)*(2*m32*mQ32*Xt*(-6*mQ32 + pow2(Xt)) - 16*pow3(
        Xt)*pow4(mQ3) + pow3(m3)*(-6*mQ32*pow2(Xt) + 6*pow4(mQ3) + pow4(Xt)) +
        mQ32*pow5(Xt) - 12*Xt*pow6(mQ3) + m3*(18*pow2(Xt)*pow4(mQ3) - 2*mQ32*
        pow4(Xt) + 3*pow6(mQ3))) + (m32 - mQ32)*(-8*(-1 + lmQ3MR)*mQ32*Xt*(6*
        mQ32 - pow2(Xt))*pow3(m3) + 4*(-1 + lmQ3MR)*pow4(m3)*(-6*mQ32*pow2(Xt)
        + 6*pow4(mQ3) + pow4(Xt)) - pow4(mQ3)*(12*(-1 + 6*lmQ3MR)*mQ32*pow2(Xt)
        + 6*(3 - 4*lmQ3MR + 6*pow2(lmQ3MR))*pow4(mQ3) + (1 - 6*lmQ3MR)*pow4(Xt)
        ) + m3*(-8*(-8 + lmQ3MR)*pow3(Xt)*pow4(mQ3) - 4*mQ32*pow5(Xt) + 48*(-1
        + lmQ3MR)*Xt*pow6(mQ3)) + m32*(12*(-7 + 8*lmQ3MR)*pow2(Xt)*pow4(mQ3) +
        (9 - 10*lmQ3MR)*mQ32*pow4(Xt) + 6*(9 - 8*lmQ3MR + 6*pow2(lmQ3MR))*pow6(
        mQ3)))))/(3.*pow2(m32 - mQ32)*pow6(mQ3));
      }
      case Limits::MQ3_EQ_M3:{
         return 4*(-23 + (8*(3*mQ32 - 3*mU32 - 8*mQ3*Xt))/(mQ32 - mU32) - 24*pow2(lmQ3MR)
        - (16*pow2(Xt))/(mQ32 - mU32) - (16*pow2(Xt))/pow2(mU3) - (128*mQ3*
        pow3(Xt))/pow2(mQ32 - mU32) + (-19*mQ32*mU32 + 8*pow4(mQ3) + 7*pow4(
        mU3))/(mQ32*mU32 - pow4(mU3)) - (10*pow4(Xt))/pow2(mQ32 - mU32) + (2*(
        27*mQ32*mU32 + 4*pow4(mQ3) - 27*pow4(mU3))*pow4(Xt))/(pow2(mU3)*pow3(
        mQ32 - mU32)) - (8*dilog(1 - mQ32/pow2(mU3))*(8*Xt*pow3(mQ3) + 4*mQ3*
        pow3(Xt) + 2*pow4(mQ3) + pow4(Xt)))/pow2(mQ32 - mU32) + (log(mU32/pow2(
        mQ3))*(-2*pow2(mQ32 - mU32)*pow2(Xt) + pow3(mQ32 - mU32) + (mQ32 +
        mU32)*pow4(Xt)))/pow3(mQ32 - mU32) + 8*lmQ3MR*(1 - mQ32/pow2(mU3) + (2*
        pow2(Xt))/pow2(mU3) + (8*mQ3*pow3(Xt))/pow2(mQ32 - mU32) - ((mQ32 + 7*
        mU32)*pow4(Xt))/pow2(-(mQ32*mU3) + pow3(mU3)) - (log(mU32/pow2(mQ3))*(
        4*mQ3*Xt*pow2(mQ32 - mU32) - 6*pow2(mQ32 - mU32)*pow2(Xt) + 3*pow3(mQ32
        - mU32) - 4*mQ3*mU32*pow3(Xt) - 4*pow3(mQ3)*pow3(Xt) + 5*mQ32*pow4(Xt)
        + 3*mU32*pow4(Xt)))/pow3(mQ32 - mU32)) + (32*mQ3*pow5(Xt))/pow3(mQ32 -
        mU32) + (log(mU32/pow2(mQ3))*(-6*(9*mQ32 - 5*mU32)*pow2(mQ32 - mU32)*
        pow2(Xt) + 32*mQ3*Xt*pow3(mQ32 - mU32) - 32*mQ3*(mQ32 - mU32)*(mQ32 +
        3*mU32)*pow3(Xt) + pow2(mQ32 - mU32)*(-26*mQ32*mU32 + 27*pow4(mQ3) +
        11*pow4(mU3)) + (36*mQ32*mU32 + 43*pow4(mQ3) - 63*pow4(mU3))*pow4(Xt) +
        32*mQ3*mU32*pow5(Xt)))/pow4(mQ32 - mU32) - (2*log(mU32/pow2(mQ3))*(-8*
        mQ3*Xt*pow3(mQ32 - mU32) - 6*pow2(Xt)*pow3(mQ32 - mU32) + 16*mQ3*pow2(
        mQ32 - mU32)*pow3(Xt) + 3*pow4(mQ32 - mU32) + 3*(pow4(mQ3) - pow4(mU3))
        *pow4(Xt) - 8*mQ3*(mQ32 + mU32)*pow5(Xt)))/pow4(mQ32 - mU32) - (4*pow2(
        log(mU32/pow2(mQ3)))*(4*Xt*(-(mQ3*mU32) + 2*pow3(mQ3))*pow3(mQ32 -
        mU32) + 4*mQ3*(mQ32 + mU32)*pow2(mQ32 - mU32)*pow3(Xt) + pow3(mQ32 -
        mU32)*(-4*mQ32*mU32 + 3*pow4(mQ3) + 2*pow4(mU3)) - 2*pow2(mQ32 - mU32)*
        pow2(Xt)*(-6*mQ32*mU32 + pow4(mQ3) + 4*pow4(mU3)) - 4*mQ3*mU32*(mQ32 +
        mU32)*pow5(Xt) + pow4(Xt)*(-10*mU32*pow4(mQ3) + mQ32*pow4(mU3) + 2*
        pow6(mQ3) + 5*pow6(mU3))))/pow5(mQ32 - mU32));
      }
      case Limits::MU3_EQ_M3:{
         return (-8*(4*mQ32*dilog(1 - mU32/pow2(mQ3))*pow3(mQ32 - mU32)*(8*Xt*pow3(mU3) +
        4*mU3*pow3(Xt) + 2*pow4(mU3) + pow4(Xt)) + pow2(mQ32 - mU32)*(-4*(-1 +
        lmQ3MR)*pow2(mU32 - pow2(Xt))*pow4(mU3) + pow4(mQ3)*(8*(-1 + 3*lmQ3MR)*
        mU32*pow2(Xt) + 64*Xt*pow3(mU3) - 32*(-2 + lmQ3MR)*mU3*pow3(Xt) + (25 -
        24*lmQ3MR + 36*pow2(lmQ3MR))*pow4(mU3) + 2*(-11 + 14*lmQ3MR)*pow4(Xt))
        + mQ32*mU3*(-8*(-2 + 3*lmQ3MR)*pow2(Xt)*pow3(mU3) + 32*(-2 + lmQ3MR)*
        mU32*pow3(Xt) - 32*Xt*pow4(mU3) + 2*(11 - 12*lmQ3MR)*mU3*pow4(Xt) + (-
        17 + 16*lmQ3MR - 12*pow2(lmQ3MR))*pow5(mU3) + 16*pow5(Xt)) - (32*mU3*Xt
        + mU32*(15 - 16*lmQ3MR + 36*pow2(lmQ3MR)) + 8*lmQ3MR*pow2(Xt))*pow6(
        mQ3) + (3 - 4*lmQ3MR + 12*pow2(lmQ3MR))*pow8(mQ3)) + 2*mQ32*pow2(log(
        mU32/pow2(mQ3)))*(mU3*pow4(mQ3)*(-38*pow2(Xt)*pow3(mU3) + 12*mU32*pow3(
        Xt) + 12*Xt*pow4(mU3) - 3*mU3*pow4(Xt) - 23*pow5(mU3) + 4*pow5(Xt)) +
        mQ32*pow3(mU3)*(32*pow2(Xt)*pow3(mU3) + 12*mU32*pow3(Xt) - 4*Xt*pow4(
        mU3) - 4*mU3*pow4(Xt) + 13*pow5(mU3) + 4*pow5(Xt)) + (20*mU32*pow2(Xt)
        - 12*Xt*pow3(mU3) - 12*mU3*pow3(Xt) + 21*pow4(mU3) + pow4(Xt))*pow6(
        mQ3) - (10*mU32*pow2(Xt) + 12*mU3*pow3(Xt) + 3*pow4(mU3) - 8*pow4(Xt))*
        pow6(mU3) - 2*(5*mU32 - 2*mU3*Xt + 2*pow2(Xt))*pow8(mQ3) + 2*power10(
        mQ3)) + (mQ32 - mU32)*log(mU32/pow2(mQ3))*(mQ32*pow3(mU3)*(2*(5 + 12*
        lmQ3MR)*pow2(Xt)*pow3(mU3) + 16*(-4 + lmQ3MR)*mU32*pow3(Xt) - 8*(-3 +
        2*lmQ3MR)*Xt*pow4(mU3) + (39 - 20*lmQ3MR)*mU3*pow4(Xt) + 3*(-3 + 4*
        lmQ3MR)*pow5(mU3) + 8*pow5(Xt)) + 2*mU3*pow4(mQ3)*(3*(1 - 12*lmQ3MR)*
        pow2(Xt)*pow3(mU3) + 32*mU32*pow3(Xt) + 12*(-3 + 2*lmQ3MR)*Xt*pow4(mU3)
        + (-17 + 4*lmQ3MR)*mU3*pow4(Xt) + (5 - 24*lmQ3MR)*pow5(mU3) + 12*pow5(
        Xt)) + (2*(-5 + 36*lmQ3MR)*mU32*pow2(Xt) + 24*(3 - 2*lmQ3MR)*Xt*pow3(
        mU3) - 16*lmQ3MR*mU3*pow3(Xt) + 2*(-5 + 36*lmQ3MR)*pow4(mU3) + (-1 +
        12*lmQ3MR)*pow4(Xt))*pow6(mQ3) + 4*pow2(mU32 - pow2(Xt))*pow6(mU3) + 2*
        ((3 - 24*lmQ3MR)*mU32 + 4*(-3 + 2*lmQ3MR)*mU3*Xt + (1 - 12*lmQ3MR)*
        pow2(Xt))*pow8(mQ3) + (-1 + 12*lmQ3MR)*power10(mQ3))))/(pow2(mQ3)*pow5(
        mQ32 - mU32));
      }
      case Limits::MQ3_EQ_MU3_EQ_M3:{
         return (-8*(-2*lmQ3MR*mQ3*Xt*(-24*mQ32*Xt - 4*mQ3*pow2(Xt) + 24*pow3(mQ3) +
        pow3(Xt)) + Xt*(-28*mQ32*pow2(Xt) + 12*Xt*pow3(mQ3) - mQ3*pow3(Xt) +
        24*pow4(mQ3) + 2*pow4(Xt)) + 36*pow2(lmQ3MR)*pow5(mQ3)))/(3.*pow5(mQ3));
      }
      default:
         break;
   };

   throw std::runtime_error("Mass limit not included!");
}

/**
 * Returns delta lambda_atas2 in the MSbar scheme for a given mass limit
 * @param limit an integer key for a mass limit
 * @param omitLogs an integer key to omit all mu terms
 * @return delta lambda_atas2 in the MSbar scheme for a given mass limit
 */
double ThresholdCalculator::getDeltaLambdaAlphatAlphas2(Limits limit, int omitLogs) const
{
   using himalaya::dilog;

   const double Mst12 = pow2(p.MSt(0));
   const double m32 = pow2(p.MG);
   const double MR2 = pow2(p.scale);
   const double mU32 = p.mu2(2,2);
   const double mQ32 = p.mq2(2,2);
   const double m3 = std::sqrt(m32);
   const double mU3 = std::sqrt(mU32);
   const double mQ3 = std::sqrt(mQ32);
   const double Xt = p.Au(2,2) - p.mu * p.vd / p.vu;
   const double lmQ3MR = omitLogs * std::log(Mst12 / MR2) + std::log(mQ32 / Mst12);
   const double dlatas2Const = 0.;
   const double z3 = 1.202056903159594;
   const double Xt4 = xtOrderLambdaAtAs2 < 4 ? 0 : pow4(Xt);
   const double Xt5 = xtOrderLambdaAtAs2 < 5 ? 0 : pow5(Xt);

   if (limit == Limits::DEGENERATE)
      limit = Limits::MQ3_EQ_MU3_EQ_M3;

   switch (limit) {
      case Limits::GENERAL: {
         const double logmqQ3 = std::log(msq2/mQ32);
         const double logU3Q3 = std::log(mU32/mQ32);
         const double logm3Q3 = std::log(m32/mQ32);

      return (3*dlatas2Const - 384*lmQ3MR*z3 + (160*lmQ3MR*m32*logmqQ3)/mQ32 +
        (160*lmQ3MR*m32*logmqQ3)/mU32 - 6096*pow2(lmQ3MR) + (2352*m32*
        pow2(lmQ3MR))/mQ32 + (2352*m32*pow2(lmQ3MR))/mU32 + (16*(-64*m32 + 207*
        (mQ32 - mU32))*logU3Q3*pow2(lmQ3MR))/(mQ32 - mU32) + (160*
        lmQ3MR*m32*(mQ32 + mU32)*Xt4*logmqQ3)/(mQ32*mU32*pow2(mQ32 -
        mU32)) + 2208*pow3(lmQ3MR) - (256*lmQ3MR*Xt5*pow3(m3))/((m32 - mQ32)*(
        m32 - mU32)*pow2(mQ32 - mU32)) + (32*lmQ3MR*(-142*m32 + 69*(mQ32 +
        mU32))*Xt4*dilog(1 - m32/mQ32))/pow3(mQ32 - mU32) + (32*lmQ3MR*(142*m32
        - 69*(mQ32 + mU32))*Xt4*dilog(1 - m32/mU32))/pow3(mQ32 - mU32) + (320*
        lmQ3MR*m32*Xt4*logmqQ3*logU3Q3)/pow3(mQ32 - mU32) + (48*(
        98*m32 + 69*(mQ32 + mU32))*Xt4*logU3Q3*pow2(lmQ3MR))/pow3(mQ32 -
        mU32) + (64*lmQ3MR*dilog(1 - m32/mQ32)*pow4(m3)*(105*mQ32*mU32 - m32*(
        27*mQ32 + 101*mU32) + 64*pow4(m3) - 41*pow4(mQ3)))/((mQ32 - mU32)*pow3(
        m32 - mQ32)) - (256*pow2(lmQ3MR)*pow4(m3))/pow4(mQ3) + (64*lmQ3MR*
        dilog(1 - m32/mU32)*pow4(m3)*(105*mQ32*mU32 - m32*(101*mQ32 + 27*mU32)
        + 64*pow4(m3) - 41*pow4(mU3)))/((-mQ32 + mU32)*pow3(m32 - mU32)) - (
        256*pow2(lmQ3MR)*pow4(m3))/pow4(mU3) + (16*Xt4*pow2(lmQ3MR)*(147*m32*
        mQ32*mU32*(mQ32 + mU32) - 16*pow2(mQ32 - mU32)*pow4(m3) + 414*pow4(mQ3)
        *pow4(mU3)))/(pow2(mQ32 - mU32)*pow4(mQ3)*pow4(mU3)) - (128*lmQ3MR*Xt5*
        pow2(logU3Q3)*(-39*m3*(mQ32 + mU32)*pow4(mU3) + 2*pow3(m3)*(15*
        mQ32*mU32 + 7*pow4(mU3)) + 8*(mQ32 + 3*mU32)*pow5(m3)))/(pow2(m32 -
        mU32)*pow4(mQ32 - mU32)) - (128*lmQ3MR*Xt5*logm3Q3*logU3Q3
        *pow3(m3)*(48*mQ32*mU32*(mQ32 + mU32) + 30*(mQ32 + mU32)*pow4(m3) -
        m32*(78*mQ32*mU32 + 47*pow4(mQ3) + 47*pow4(mU3)) + 16*pow6(m3)))/(pow2(
        m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) - (256*lmQ3MR*Xt5*
        logm3Q3*pow3(m3)*(-43*m32*mQ32*mU32*(mQ32 + mU32) + pow4(m3)*(38*
        mQ32*mU32 - 4*pow4(mQ3) - 4*pow4(mU3)) + 48*pow4(mQ3)*pow4(mU3) + 4*(
        mQ32 + mU32)*pow6(m3)))/(mQ32*mU32*pow2(m32 - mQ32)*pow2(m32 - mU32)*
        pow2(mQ32 - mU32)) + (128*lmQ3MR*Xt5*logU3Q3*(78*m3*pow4(mQ3)*
        pow4(mU3) - pow3(m3)*(53*mU32*pow4(mQ3) + 71*mQ32*pow4(mU3)) + (45*
        mQ32*mU32 - 23*pow4(mQ3) - 8*pow4(mU3))*pow5(m3) + 8*(3*mQ32 + mU32)*
        pow7(m3)))/(mQ32*(-m32 + mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) + (
        16*lmQ3MR*pow2(logU3Q3)*(pow4(m3)*(519*mQ32*mU32 - 263*pow4(mU3)
        ) - (239*mQ32 + 81*mU32)*pow6(m3) + 138*(mQ32 - mU32)*pow6(mU3) + m32*(
        -414*mQ32*pow4(mU3) + 350*pow6(mU3)) + 128*pow8(m3)))/((-mQ32 + mU32)*
        pow3(m32 - mU32)) + 64*lmQ3MR*Xt*((138*lmQ3MR*m3*logU3Q3)/(mQ32
        - mU32) + (10*m3*logmqQ3*logU3Q3)/(mQ32 - mU32) - (64*
        pow3(m3))/(mQ32*mU32) + (32*lmQ3MR*pow3(m3))/(mQ32*mU32) - (8*dilog(1 -
        m32/mQ32)*(-48*mQ32*pow3(m3) + 47*pow5(m3)))/((mQ32 - mU32)*pow2(m32 -
        mQ32)) - (8*dilog(1 - m32/mU32)*(-48*mU32*pow3(m3) + 47*pow5(m3)))/((-
        mQ32 + mU32)*pow2(m32 - mU32)) + (8*logm3Q3*pow5(m3)*(-11*mQ32*
        mU32*(mQ32 + mU32) - 18*(mQ32 + mU32)*pow4(m3) + 10*m32*(3*mQ32*mU32 +
        pow4(mQ3) + pow4(mU3)) + 8*pow6(m3)))/(mQ32*mU32*pow2(m32 - mQ32)*pow2(
        m32 - mU32)) + (m3*pow2(logU3Q3)*((189*mQ32 - 221*mU32)*pow4(m3)
        + 111*mQ32*pow4(mU3) + m32*(-318*mQ32*mU32 + 382*pow4(mU3)) - 143*pow6(
        mU3)))/(pow2(m32 - mU32)*pow2(mQ32 - mU32)) + (2*logU3Q3*(-107*
        m3*pow4(mQ3)*pow4(mU3) + pow3(m3)*(263*mU32*pow4(mQ3) + 124*mQ32*pow4(
        mU3)) - (289*mQ32*mU32 + 152*pow4(mQ3) + 8*pow4(mU3))*pow5(m3) + (161*
        mQ32 + 8*mU32)*pow7(m3)))/(mQ32*(-m32 + mQ32)*(mQ32 - mU32)*pow2(m32 -
        mU32)) - (6*m3*logm3Q3*logU3Q3*(-12*m32*mQ32*mU32*(mQ32 +
        mU32) - 2*pow4(mQ3)*pow4(mU3) + pow4(m3)*(56*mQ32*mU32 + 11*pow4(mQ3) +
        11*pow4(mU3)) - 38*(mQ32 + mU32)*pow6(m3) + 24*pow8(m3)))/((mQ32 -
        mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32))) - (16*lmQ3MR*Xt4*pow2(logU3Q3)*
        (4*(89*mQ32*mU32 + 8*pow4(mQ3) - 311*pow4(mU3))*pow6(m3) +
        69*(-4*mQ32*mU32 + pow4(mQ3) - 5*pow4(mU3))*pow6(mU3) + 2*pow4(m3)*(51*
        mU32*pow4(mQ3) - 310*mQ32*pow4(mU3) + 455*pow6(mU3)) - 16*(7*mQ32 - 23*
        mU32)*pow8(m3) + m32*(-207*pow4(mQ3)*pow4(mU3) + 652*mQ32*pow6(mU3) +
        315*pow8(mU3))))/(pow3(m32 - mU32)*pow4(mQ32 - mU32)) + (16*lmQ3MR*
        logm3Q3*logU3Q3*pow4(m3)*((372*mQ32*mU32 + 454*pow4(mQ3) +
        454*pow4(mU3))*pow6(m3) + 41*pow4(mU3)*pow6(mQ3) + 41*pow4(mQ3)*pow6(
        mU3) - pow4(m3)*(273*mU32*pow4(mQ3) + 273*mQ32*pow4(mU3) + 367*pow6(
        mQ3) + 367*pow6(mU3)) - 320*(mQ32 + mU32)*pow8(m3) - 105*mU32*pow8(mQ3)
        - 105*mQ32*pow8(mU3) + m32*(-246*pow4(mQ3)*pow4(mU3) + 342*mU32*pow6(
        mQ3) + 342*mQ32*pow6(mU3) + 101*pow8(mQ3) + 101*pow8(mU3)) + 128*
        power10(m3)))/((mQ32 - mU32)*pow3(m32 - mQ32)*pow3(m32 - mU32)) + (16*
        lmQ3MR*m32*Xt4*logm3Q3*logU3Q3*(88*pow12(m3) + 24*pow6(
        mQ3)*pow6(mU3) + pow6(m3)*(2221*mU32*pow4(mQ3) + 2221*mQ32*pow4(mU3) +
        567*pow6(mQ3) + 567*pow6(mU3)) - 8*(211*mQ32*mU32 + 110*pow4(mQ3) +
        110*pow4(mU3))*pow8(m3) - pow4(m3)*(1694*pow4(mQ3)*pow4(mU3) + 952*
        mU32*pow6(mQ3) + 952*mQ32*pow6(mU3) + 101*pow8(mQ3) + 101*pow8(mU3)) +
        m32*(353*pow4(mU3)*pow6(mQ3) + 353*pow4(mQ3)*pow6(mU3) + 105*mU32*pow8(
        mQ3) + 105*mQ32*pow8(mU3)) + 322*(mQ32 + mU32)*power10(m3)))/(pow3(m32
        - mQ32)*pow3(m32 - mU32)*pow3(mQ32 - mU32)) + (16*lmQ3MR*Xt4*(32*pow12(
        m3)*pow2(mQ32 - mU32) - 6*m32*(mQ32 + mU32)*pow4(mQ3)*pow4(mU3)*(mQ32*(
        10*msq2 - 383*mU32) + 10*msq2*mU32 + pow4(mQ3) + pow4(mU3)) - 2*mQ32*
        mU32*(mQ32 + mU32)*(mQ32*(30*msq2 - 1007*mU32) + 30*msq2*mU32 + 167*
        pow4(mQ3) + 167*pow4(mU3))*pow6(m3) + 3*(10*mQ32*(msq2 - 39*mU32) + 10*
        msq2*mU32 + pow4(mQ3) + pow4(mU3))*pow6(mQ3)*pow6(mU3) + mQ32*mU32*
        pow4(m3)*(2*pow4(mQ3)*(75*msq2*mU32 - 2269*pow4(mU3)) + (30*msq2 - 812*
        mU32)*pow6(mQ3) + 2*mQ32*(75*msq2*pow4(mU3) - 406*pow6(mU3)) + 3*(10*
        msq2 + mU32)*pow6(mU3) + 3*pow8(mQ3)) + pow8(m3)*(10*pow4(mQ3)*(3*msq2*
        mU32 - 67*pow4(mU3)) + 611*mU32*pow6(mQ3) + mQ32*(30*msq2*pow4(mU3) +
        611*pow6(mU3)) + 32*pow8(mQ3) + 32*pow8(mU3)) - 8*(27*mU32*pow4(mQ3) +
        27*mQ32*pow4(mU3) + 8*pow6(mQ3) + 8*pow6(mU3))*power10(m3)))/(pow2(m32
        - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32)*pow4(mQ3)*pow4(mU3)) + (8*
        lmQ3MR*(64*pow12(m3)*(pow4(mQ3) + pow4(mU3)) - 12*m32*(mQ32 + mU32)*
        pow4(mQ3)*pow4(mU3)*(10*msq2*mU32 + mQ32*(10*msq2 + 147*mU32) + pow4(
        mQ3) + pow4(mU3)) - 4*mQ32*mU32*(mQ32 + mU32)*(30*msq2*mU32 + 6*mQ32*(
        5*msq2 + 154*mU32) + 167*pow4(mQ3) + 167*pow4(mU3))*pow6(m3) + (6*mU32*
        (10*msq2 + mU32) + mQ32*(60*msq2 + 559*mU32) + 6*pow4(mQ3))*pow6(mQ3)*
        pow6(mU3) + mQ32*mU32*pow4(m3)*(12*pow4(mQ3)*(25*msq2*mU32 + 408*pow4(
        mU3)) + (60*msq2 + 1829*mU32)*pow6(mQ3) + 6*(10*msq2 + mU32)*pow6(mU3)
        + mQ32*(300*msq2*pow4(mU3) + 1829*pow6(mU3)) + 6*pow8(mQ3)) + pow8(m3)*
        (pow4(mQ3)*(60*msq2*mU32 + 3179*pow4(mU3)) + 1350*mU32*pow6(mQ3) + 30*
        mQ32*(2*msq2*pow4(mU3) + 45*pow6(mU3)) + 64*pow8(mQ3) + 64*pow8(mU3)) -
        16*(43*mU32*pow4(mQ3) + 43*mQ32*pow4(mU3) + 8*pow6(mQ3) + 8*pow6(mU3))*
        power10(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)*pow4(mQ3)*pow4(mU3)) +
        (16*lmQ3MR*Xt4*logU3Q3*(pow12(m3)*(-784*mQ32*mU32 + 33*pow4(mQ3)
        - 65*pow4(mU3)) - 2*m32*pow4(mQ3)*(-60*msq2*mU32 + mQ32*(-90*msq2 +
        1773*mU32) + 744*pow4(mQ3) + 1283*pow4(mU3))*pow6(mU3) + mQ32*pow4(m3)*
        pow4(mU3)*((-180*msq2 + 3613*mU32)*pow4(mQ3) - 60*msq2*pow4(mU3) +
        mQ32*(-360*msq2*mU32 + 5687*pow4(mU3)) + 1899*pow6(mQ3) + 1209*pow6(
        mU3)) + pow8(m3)*(-60*pow4(mQ3)*(2*msq2*mU32 + 21*pow4(mU3)) + 1111*
        mU32*pow6(mQ3) + mQ32*(-180*msq2*pow4(mU3) + 75*pow6(mU3)) + 33*pow8(
        mQ3) - 243*pow8(mU3)) + 12*(37*mQ32 - 5*msq2 + 106*mU32)*pow6(mQ3)*
        pow8(mU3) + mU32*pow6(m3)*(3*pow4(mQ3)*(120*msq2*mU32 - 809*pow4(mU3))
        + 4*(15*msq2 - 602*mU32)*pow6(mQ3) + 6*mQ32*(30*msq2*pow4(mU3) - 389*
        pow6(mU3)) - 880*pow8(mQ3) + 81*pow8(mU3)) + (557*mU32*pow4(mQ3) + 2*
        mQ32*(30*msq2*mU32 + 913*pow4(mU3)) - 66*pow6(mQ3) + 227*pow6(mU3))*
        power10(m3)))/(mQ32*mU32*pow2(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 -
        mU32)) + (64*lmQ3MR*pow3(Xt)*(-10*m3*(mQ32 - mU32)*(mQ32 + mU32)*logmqQ3*
        logU3Q3 - 20*m3*logmqQ3*pow2(mQ32 - mU32) + 4*
        (mQ32 - mU32)*dilog(1 - m32/mQ32)*(31*m3*(mQ32 + mU32) - 64*pow3(m3)) -
        2*lmQ3MR*(mQ32 - mU32)*logU3Q3*(69*m3*(mQ32 + mU32) + 32*pow3(
        m3)) + 4*(mQ32 - mU32)*dilog(1 - m32/mU32)*(-31*m3*(mQ32 + mU32) + 64*
        pow3(m3)) - (4*lmQ3MR*pow2(mQ32 - mU32)*(69*m3*mQ32*mU32 + 8*(mQ32 +
        mU32)*pow3(m3)))/(mQ32*mU32) + (4*pow2(mQ32 - mU32)*(193*mQ32*mU32*(
        mQ32 + mU32)*pow3(m3) - 201*m3*pow4(mQ3)*pow4(mU3) + (-185*mQ32*mU32 +
        16*pow4(mQ3) + 16*pow4(mU3))*pow5(m3) - 16*(mQ32 + mU32)*pow7(m3)))/(
        mQ32*(-m32 + mQ32)*(m32 - mU32)*mU32) - (4*m3*(mQ32 - mU32)*logU3Q3*(
        (101*mQ32*mU32 + 44*pow4(mQ3) + 8*pow4(mU3))*pow6(m3) - 53*pow4(
        mU3)*pow6(mQ3) - 180*pow4(mQ3)*pow6(mU3) - pow4(m3)*(224*mU32*pow4(mQ3)
        + 332*mQ32*pow4(mU3) + 71*pow6(mQ3) + 8*pow6(mU3)) + m32*(376*pow4(mQ3)
        *pow4(mU3) + 116*mU32*pow6(mQ3) + 191*mQ32*pow6(mU3)) + 32*mQ32*pow8(
        m3)))/(mQ32*(-m32 + mQ32)*pow2(m32 - mU32)) + (m3*pow2(logU3Q3)*
        (-251*pow4(mQ3)*pow4(mU3) + pow4(m3)*(-256*mQ32*mU32 - 47*pow4(mQ3) +
        431*pow4(mU3)) + 128*(mQ32 - mU32)*pow6(m3) + m32*(310*mU32*pow4(mQ3) -
        240*mQ32*pow4(mU3) - 326*pow6(mU3)) + 376*mQ32*pow6(mU3) + 3*pow8(mU3))
        )/pow2(m32 - mU32) - (4*m3*(mQ32 - mU32)*logm3Q3*logU3Q3*(
        3*(mQ32 + mU32)*pow4(mQ3)*pow4(mU3) + 4*(112*mQ32*mU32 + 9*pow4(mQ3) +
        9*pow4(mU3))*pow6(m3) + pow4(m3)*(-220*mU32*pow4(mQ3) - 220*mQ32*pow4(
        mU3) + 54*pow6(mQ3) + 54*pow6(mU3)) - 18*m32*(-10*pow4(mQ3)*pow4(mU3) +
        3*mU32*pow6(mQ3) + 3*mQ32*pow6(mU3)) - 181*(mQ32 + mU32)*pow8(m3) + 96*
        power10(m3)))/(pow2(m32 - mQ32)*pow2(m32 - mU32)) - (4*m3*logm3Q3
        *pow2(mQ32 - mU32)*(5*m32*(mQ32 + mU32)*pow4(mQ3)*pow4(mU3) + 6*pow6(
        mQ3)*pow6(mU3) + pow6(m3)*(47*mU32*pow4(mQ3) + 47*mQ32*pow4(mU3) + 24*
        pow6(mQ3) + 24*pow6(mU3)) - pow4(m3)*(28*pow4(mQ3)*pow4(mU3) + 27*mU32*
        pow6(mQ3) + 27*mQ32*pow6(mU3)) - 4*(7*mQ32*mU32 + 10*pow4(mQ3) + 10*
        pow4(mU3))*pow8(m3) + 16*(mQ32 + mU32)*power10(m3)))/(mQ32*mU32*pow2(
        m32 - mQ32)*pow2(m32 - mU32))))/pow4(mQ32 - mU32) + (32*lmQ3MR*m32*Xt4*
        logm3Q3*(16*pow14(m3)*pow2(mQ32 - mU32) + pow4(m3)*pow4(mQ3)*
        pow4(mU3)*(89*mU32*pow4(mQ3) + 89*mQ32*pow4(mU3) - 184*pow6(mQ3) - 184*
        pow6(mU3)) + 2*m32*(49*mQ32*mU32 + 69*pow4(mQ3) + 69*pow4(mU3))*pow6(
        mQ3)*pow6(mU3) - 3*pow12(m3)*(15*mU32*pow4(mQ3) + 15*mQ32*pow4(mU3) +
        16*pow6(mQ3) + 16*pow6(mU3)) - 75*(mQ32 + mU32)*pow8(mQ3)*pow8(mU3) +
        mQ32*mU32*pow6(m3)*(-692*pow4(mQ3)*pow4(mU3) + 219*mU32*pow6(mQ3) +
        219*mQ32*pow6(mU3) + 133*pow8(mQ3) + 133*pow8(mU3)) + (-192*pow4(mQ3)*
        pow4(mU3) + 311*mU32*pow6(mQ3) + 311*mQ32*pow6(mU3) + 48*pow8(mQ3) +
        48*pow8(mU3))*power10(m3) - pow8(m3)*(-190*pow4(mU3)*pow6(mQ3) - 190*
        pow4(mQ3)*pow6(mU3) + 367*mU32*pow8(mQ3) + 367*mQ32*pow8(mU3) + 16*
        power10(mQ3) + 16*power10(mU3))))/(pow2(mQ32 - mU32)*pow3(-m32 + mQ32)*
        pow3(m32 - mU32)*pow4(mQ3)*pow4(mU3)) - (16*lmQ3MR*m32*logm3Q3*(-
        32*pow14(m3)*(pow4(mQ3) + pow4(mU3)) + 3*m32*(-24*mQ32*mU32 + 23*pow4(
        mQ3) + 23*pow4(mU3))*pow6(mQ3)*pow6(mU3) + pow12(m3)*(346*mU32*pow4(
        mQ3) + 346*mQ32*pow4(mU3) + 96*pow6(mQ3) + 96*pow6(mU3)) + pow4(m3)*
        pow4(mQ3)*pow4(mU3)*(-155*mU32*pow4(mQ3) - 155*mQ32*pow4(mU3) + 237*
        pow6(mQ3) + 237*pow6(mU3)) + 12*(mQ32 + mU32)*pow8(mQ3)*pow8(mU3) -
        mQ32*mU32*pow6(m3)*(-318*pow4(mQ3)*pow4(mU3) + 739*mU32*pow6(mQ3) +
        739*mQ32*pow6(mU3) + 330*pow8(mQ3) + 330*pow8(mU3)) - 2*(285*pow4(mQ3)*
        pow4(mU3) + 503*mU32*pow6(mQ3) + 503*mQ32*pow6(mU3) + 48*pow8(mQ3) +
        48*pow8(mU3))*power10(m3) + pow8(m3)*(738*pow4(mU3)*pow6(mQ3) + 738*
        pow4(mQ3)*pow6(mU3) + 990*mU32*pow8(mQ3) + 990*mQ32*pow8(mU3) + 32*
        power10(mQ3) + 32*power10(mU3))))/(pow3(-m32 + mQ32)*pow3(m32 - mU32)*
        pow4(mQ3)*pow4(mU3)) + (16*lmQ3MR*logU3Q3*(pow12(m3)*(96*mQ32*
        mU32 + 33*pow4(mQ3) - pow4(mU3)) + mQ32*pow4(m3)*pow4(mU3)*(1491*mU32*
        pow4(mQ3) - 2180*mQ32*pow4(mU3) + 1843*pow6(mQ3) - 514*pow6(mU3)) +
        pow8(m3)*(990*pow4(mQ3)*pow4(mU3) + 2083*mU32*pow6(mQ3) - 1775*mQ32*
        pow6(mU3) + 33*pow8(mQ3) - 51*pow8(mU3)) + 381*(mQ32 - mU32)*pow6(mQ3)*
        pow8(mU3) + (-1233*mU32*pow4(mQ3) + 624*mQ32*pow4(mU3) - 66*pow6(mQ3) +
        35*pow6(mU3))*power10(m3) + pow6(m3)*(-3472*pow4(mU3)*pow6(mQ3) + 1568*
        pow4(mQ3)*pow6(mU3) - 942*mU32*pow8(mQ3) + 1549*mQ32*pow8(mU3) + 17*
        power10(mU3)) + m32*(-1295*pow6(mU3)*pow8(mQ3) + 285*pow6(mQ3)*pow8(
        mU3) + 882*pow4(mQ3)*power10(mU3))))/(mQ32*(mQ32 - mU32)*mU32*pow2(m32
        - mQ32)*pow3(m32 - mU32)) + 32*lmQ3MR*pow2(Xt)*((118*m32)/(mQ32*mU32) -
        (3*(-8*m32 + mQ32 + 10*msq2 + mU32))/(mQ32*mU32) - (10*m32*logmqQ3)/
        (mQ32*mU32) - (768*m32)/pow2(mQ32 - mU32) - (192*m32*(-2*m32 +
        mQ32 + mU32)*dilog(1 - m32/mQ32))/pow3(mQ32 - mU32) + (192*m32*(-2*m32
        + mQ32 + mU32)*dilog(1 - m32/mU32))/pow3(mQ32 - mU32) + (3*lmQ3MR*
        logU3Q3*(32*m32*(mQ32 + mU32) - 69*pow2(mQ32 - mU32)))/pow3(mQ32 -
        mU32) + (32*(2*m32*mU32 - 3*pow4(m3)))/(mQ32*mU32*(-m32 + mU32)) - (64*
        pow4(m3))/(mU32*pow4(mQ3)) + (4*(8*m32 - 9*mQ32)*(m32 - mQ32 - mU32)*
        pow4(m3))/((m32 - mQ32)*(m32 - mU32)*mU32*pow4(mQ3)) - (32*pow4(m3))/(
        mQ32*pow4(mU3)) + (-46*m32*mQ32*mU32 + 62*(mQ32 + mU32)*pow4(m3) - 78*
        pow6(m3))/(mQ32*(-m32 + mQ32)*(m32 - mU32)*mU32) + (lmQ3MR*m32*(16*m32*
        (mQ32 + mU32) - (3*(-162*pow4(mQ3)*pow4(mU3) + 49*mU32*pow6(mQ3) + 49*
        mQ32*pow6(mU3)))/pow2(mQ32 - mU32)))/(pow4(mQ3)*pow4(mU3)) + (pow2(
        logU3Q3)*((364*mQ32*mU32 + 37*pow4(mQ3) + 751*pow4(mU3))*pow6(m3) -
        138*(-3*mQ32*mU32 + pow4(mQ3) + 2*pow4(mU3))*pow6(mU3) - 3*pow4(m3)*(
        103*mU32*pow4(mQ3) - 88*mQ32*pow4(mU3) + 241*pow6(mU3)) - 128*(mQ32 +
        5*mU32)*pow8(m3) + m32*(414*pow4(mQ3)*pow4(mU3) - 922*mQ32*pow6(mU3) +
        700*pow8(mU3)) + 192*power10(m3)))/(pow3(m32 - mU32)*pow3(-mQ32 + mU32)
        ) + (2*m32*logm3Q3*(117*m32*(mQ32 + mU32)*pow4(mQ3)*pow4(mU3) -
        6*pow6(mQ3)*pow6(mU3) + 2*pow6(m3)*(141*mU32*pow4(mQ3) + 141*mQ32*pow4(
        mU3) + 8*pow6(mQ3) + 8*pow6(mU3)) - pow4(m3)*(353*pow4(mQ3)*pow4(mU3) +
        133*mU32*pow6(mQ3) + 133*mQ32*pow6(mU3)) - (173*mQ32*mU32 + 32*pow4(
        mQ3) + 32*pow4(mU3))*pow8(m3) + 16*(mQ32 + mU32)*power10(m3)))/(pow2(
        m32 - mQ32)*pow2(m32 - mU32)*pow4(mQ3)*pow4(mU3)) - (logU3Q3*(-
        450*pow2(mQ32 - mU32)*pow6(mQ3)*pow6(mU3) + 2*m32*pow4(mQ3)*pow4(mU3)*(
        -653*mU32*pow4(mQ3) - 53*mQ32*pow4(mU3) + 613*pow6(mQ3) + 477*pow6(mU3)
        ) - 2*pow8(m3)*(-98*pow4(mQ3)*pow4(mU3) + 551*mU32*pow6(mQ3) + 1099*
        mQ32*pow6(mU3) + 33*pow8(mQ3) - 49*pow8(mU3)) - mQ32*mU32*pow4(m3)*(-
        2316*pow4(mQ3)*pow4(mU3) + 1310*mU32*pow6(mQ3) + 2810*mQ32*pow6(mU3) +
        817*pow8(mQ3) + 451*pow8(mU3)) + (93*mU32*pow4(mQ3) + 675*mQ32*pow4(
        mU3) + 33*pow6(mQ3) - 33*pow6(mU3))*power10(m3) + pow6(m3)*(-803*pow4(
        mU3)*pow6(mQ3) + 1695*pow4(mQ3)*pow6(mU3) + 1830*mU32*pow8(mQ3) + 1902*
        mQ32*pow8(mU3) + 33*power10(mQ3) - 49*power10(mU3))))/(mQ32*mU32*pow2(
        m32 - mQ32)*pow2(m32 - mU32)*pow3(mQ32 - mU32)) + (logm3Q3*logU3Q3*
        pow4(m3)*(384*pow12(m3) + 24*pow6(m3)*(-177*mU32*pow4(mQ3) -
        177*mQ32*pow4(mU3) + 17*pow6(mQ3) + 17*pow6(mU3)) + (4244*mQ32*mU32 +
        758*pow4(mQ3) + 758*pow4(mU3))*pow8(m3) + pow4(m3)*(5262*pow4(mQ3)*
        pow4(mU3) + 744*mU32*pow6(mQ3) + 744*mQ32*pow6(mU3) - 495*pow8(mQ3) -
        495*pow8(mU3)) + 3*mQ32*mU32*(186*pow4(mQ3)*pow4(mU3) + 6*mU32*pow6(
        mQ3) + 6*mQ32*pow6(mU3) - 35*pow8(mQ3) - 35*pow8(mU3)) - 1152*(mQ32 +
        mU32)*power10(m3) + m32*(-1750*pow4(mU3)*pow6(mQ3) - 1750*pow4(mQ3)*
        pow6(mU3) + 497*mU32*pow8(mQ3) + 497*mQ32*pow8(mU3) + 101*power10(mQ3)
        + 101*power10(mU3))))/(pow3(m32 - mQ32)*pow3(m32 - mU32)*pow3(mQ32 -
        mU32))))/3.;
      }

      case Limits::MQ3_EQ_MU3: {
         const double logmqQ3 = std::log(msq2/mQ32);
         const double logm3Q3 = std::log(m32/mQ32);

         return -(384*lmQ3MR*dilog(1 - m32/mQ32)*pow3(m3)*(-137*m3*mQ32 - 184*m32*Xt +
        192*mQ32*Xt + 133*pow3(m3))*pow8(mQ3) + (m32 - mQ32)*(160*lmQ3MR*m3*
        mQ32*logmqQ3*pow2(m32 - mQ32)*(2*mQ32*Xt*(-6*mQ32 + pow2(Xt)) +
        m3*(Xt4 - 6*mQ32*pow2(Xt) + 6*pow4(mQ3))) + 9*dlatas2Const*pow2(m32 -
        mQ32)*pow8(mQ3) + 6624*pow2(m32 - mQ32)*pow3(lmQ3MR)*pow8(mQ3) - 8*
        pow2(lmQ3MR)*pow2(m32 - mQ32)*(96*pow4(m3)*(Xt4 - 4*mQ32*pow2(Xt) + 2*
        pow4(mQ3)) - 256*pow3(m3)*(-(mQ32*pow3(Xt)) + 3*Xt*pow4(mQ3)) - 6*m32*(
        49*mQ32*Xt4 - 326*pow2(Xt)*pow4(mQ3) + 358*pow6(mQ3)) + 552*m3*(-(pow3(
        Xt)*pow4(mQ3)) + 6*Xt*pow6(mQ3)) + 9*(23*Xt4*pow4(mQ3) - 276*pow2(Xt)*
        pow6(mQ3) + 254*pow8(mQ3))) + 8*lmQ3MR*(104*m3*(3*Xt5 - 49*mQ32*pow3(
        Xt) + 42*Xt*pow4(mQ3))*pow6(mQ3) - 16*pow6(m3)*(61*mQ32*Xt4 - 354*pow2(
        Xt)*pow4(mQ3) + 306*pow6(mQ3)) + 8*pow5(m3)*(8*mQ32*Xt5 - 249*pow3(Xt)*
        pow4(mQ3) + 1302*Xt*pow6(mQ3)) + 3*pow6(mQ3)*(20*msq2*Xt4 - 12*mQ32*(-
        9*Xt4 + 10*msq2*pow2(Xt)) + 24*(5*msq2 - 53*pow2(Xt))*pow4(mQ3) + (571
        - 48*z3)*pow6(mQ3)) + 2*m32*pow4(mQ3)*(-60*msq2*Xt4 + mQ32*(-926*Xt4 +
        360*msq2*pow2(Xt)) + (-360*msq2 + 9336*pow2(Xt))*pow4(mQ3) + 3*(-1345 +
        48*z3)*pow6(mQ3)) - 512*(-(mQ32*pow3(Xt)) + 3*Xt*pow4(mQ3))*pow7(m3) +
        192*(Xt4 - 4*mQ32*pow2(Xt) + 2*pow4(mQ3))*pow8(m3) - 8*pow3(m3)*(45*
        Xt5*pow4(mQ3) - 794*pow3(Xt)*pow6(mQ3) + 1680*Xt*pow8(mQ3)) + pow4(m3)*
        (60*mQ32*msq2*Xt4 - 8*(-287*Xt4 + 45*msq2*pow2(Xt))*pow4(mQ3) + 120*(3*
        msq2 - 163*pow2(Xt))*pow6(mQ3) - 9*(-1213 + 16*z3)*pow8(mQ3)))) - 32*
        lmQ3MR*m3*logm3Q3*(pow3(m3)*(-308*Xt4 + 2766*mQ32*pow2(Xt) + 399*
        pow4(mQ3))*pow6(mQ3) - 12*m32*(-8*Xt5 + 131*mQ32*pow3(Xt) + 78*Xt*pow4(
        mQ3))*pow6(mQ3) + 8*pow6(m3)*(4*mQ32*Xt5 - 93*pow3(Xt)*pow4(mQ3) + 342*
        Xt*pow6(mQ3)) - 3*(79*mQ32*Xt4 - 442*pow2(Xt)*pow4(mQ3) + 410*pow6(mQ3)
        )*pow7(m3) - 128*(-(mQ32*pow3(Xt)) + 3*Xt*pow4(mQ3))*pow8(m3) + 6*m3*(
        Xt4 - 6*mQ32*pow2(Xt) + 6*pow4(mQ3))*pow8(mQ3) + pow5(m3)*(487*Xt4*
        pow4(mQ3) - 3816*pow2(Xt)*pow6(mQ3) + 759*pow8(mQ3)) - 4*pow4(m3)*(31*
        Xt5*pow4(mQ3) - 530*pow3(Xt)*pow6(mQ3) + 372*Xt*pow8(mQ3)) + 48*(Xt4 -
        4*mQ32*pow2(Xt) + 2*pow4(mQ3))*pow9(m3) + 12*Xt*(-6*mQ32 + pow2(Xt))*
        power10(mQ3)))/(9.*pow3(-m32 + mQ32)*pow8(mQ3));
      }

      case Limits::MQ3_EQ_M3: {
         const double logmqQ3 = std::log(msq2/mQ32);
         const double logU3Q3 = std::log(mU32/mQ32);

         return (3*dlatas2Const + (8*lmQ3MR*(3857 + (216*mQ32)/mU32))/3. - (32*lmQ3MR*(
        3853*mQ32 - 2701*mU32))/(9.*(mQ32 - mU32)) - 384*lmQ3MR*z3 + 160*
        lmQ3MR*logmqQ3 + (160*lmQ3MR*mQ32*logmqQ3)/mU32 + (16*
        lmQ3MR*(569*mQ32 - 377*mU32)*logU3Q3)/(3.*(mQ32 - mU32)) - 4000*
        pow2(lmQ3MR) + (2352*mQ32*pow2(lmQ3MR))/mU32 + (16*(143*mQ32 - 207*
        mU32)*logU3Q3*pow2(lmQ3MR))/(mQ32 - mU32) + (160*lmQ3MR*(mQ32 +
        mU32)*Xt4*logmqQ3)/pow2(-(mQ32*mU3) + pow3(mU3)) + 2208*pow3(
        lmQ3MR) + (32*lmQ3MR*(73*mQ32 - 69*mU32)*Xt4*dilog(1 - mQ32/mU32))/
        pow3(mQ32 - mU32) + (320*lmQ3MR*mQ32*Xt4*logmqQ3*logU3Q3)
        /pow3(mQ32 - mU32) + (48*(167*mQ32 + 69*mU32)*Xt4*logU3Q3*pow2(
        lmQ3MR))/pow3(mQ32 - mU32) + (64*lmQ3MR*mQ3*(16*mQ32 + 201*mU32)*Xt5)/(
        mU32*pow3(-mQ32 + mU32)) + (64*lmQ3MR*(37*mQ32 - 41*mU32)*dilog(1 -
        mQ32/mU32)*pow4(mQ3))/pow3(mQ32 - mU32) + (32*lmQ3MR*Xt4*(1343*mQ32*
        mU32 + 54*pow4(mQ3) - 437*pow4(mU3)))/(3.*mU32*pow3(mQ32 - mU32)) - (
        256*pow2(lmQ3MR)*pow4(mQ3))/pow4(mU3) + (16*Xt4*pow2(lmQ3MR)*(179*mQ32*
        mU32 - 16*pow4(mQ3) + 545*pow4(mU3)))/(pow2(mQ32 - mU32)*pow4(mU3)) - (
        32*lmQ3MR*mQ3*(249*mQ32 + 185*mU32)*Xt5*logU3Q3)/pow4(mQ32 -
        mU32) - (64*lmQ3MR*Xt5*(-5*mQ3*mU32 + pow3(mQ3)))/pow4(mQ32 - mU32) + (
        16*lmQ3MR*Xt4*logU3Q3*(744*mQ32*mU32 + 1553*pow4(mQ3) - 377*
        pow4(mU3)))/(3.*pow4(mQ32 - mU32)) - (32*lmQ3MR*mQ3*Xt5*logU3Q3*
        (236*mQ32*mU32 + 97*pow4(mQ3) - 349*pow4(mU3)))/pow5(mQ32 - mU32) + 64*
        lmQ3MR*pow3(Xt)*((-20*mQ3*logmqQ3)/pow2(mQ32 - mU32) + (4*mQ3*(
        8*mQ32 - 21*mU32))/pow2(-(mQ32*mU3) + pow3(mU3)) - (4*lmQ3MR*(77*mQ3*
        mU32 + 8*pow3(mQ3)))/pow2(-(mQ32*mU3) + pow3(mU3)) + (mQ3*(157*mQ32 -
        213*mU32)*logU3Q3)/pow3(mQ32 - mU32) - (10*mQ3*(mQ32 + mU32)*
        logmqQ3*logU3Q3)/pow3(mQ32 - mU32) + (4*dilog(1 - mQ32/
        mU32)*(-31*mQ3*mU32 + 33*pow3(mQ3)))/pow3(mQ32 - mU32) - (2*lmQ3MR*
        logU3Q3*(69*mQ3*mU32 + 101*pow3(mQ3)))/pow3(mQ32 - mU32) + (4*mQ3*(
        219*mQ32*mU32 + 16*pow4(mQ3) - 227*pow4(mU3)))/(mU32*pow3(mQ32 - mU32))
        + (mQ3*logU3Q3*(330*mQ32*mU32 + 437*pow4(mQ3) - 703*pow4(mU3)))/
        pow4(mQ32 - mU32) + (mQ3*pow2(logU3Q3)*(7*mU32*pow4(mQ3) - 53*
        mQ32*pow4(mU3) + 81*pow6(mQ3) - 3*pow6(mU3)))/pow5(mQ32 - mU32)) + (16*
        lmQ3MR*pow2(logU3Q3)*(-327*mU32*pow4(mQ3) + 350*mQ32*pow4(mU3) +
        111*pow6(mQ3) - 138*pow6(mU3)))/pow3(mQ32 - mU32) + (64*lmQ3MR*mQ3*Xt*(
        2*(8*(-5 + 2*lmQ3MR)*mQ32 + (183 - 16*lmQ3MR)*mU32)*pow2(mQ32 - mU32) +
        mU32*logU3Q3*(12*(57 - 23*lmQ3MR)*mQ32*mU32 + 10*logmqQ3*
        pow2(mQ32 - mU32) + 23*(-17 + 6*lmQ3MR)*pow4(mQ3) + 3*(-95 + 46*lmQ3MR)
        *pow4(mU3)) + 8*dilog(1 - mQ32/mU32)*(47*mU32*pow4(mQ3) - 48*mQ32*pow4(
        mU3)) + pow2(logU3Q3)*(189*mU32*pow4(mQ3) - 350*mQ32*pow4(mU3) +
        143*pow6(mU3))))/(mU32*pow3(mQ32 - mU32)) - (128*lmQ3MR*Xt5*pow2(logU3Q3)*
        (-25*pow3(mQ3)*pow4(mU3) + 54*mU32*pow5(mQ3) - 39*mQ3*pow6(
        mU3) + 8*pow7(mQ3)))/pow6(mQ32 - mU32) + (8*lmQ3MR*Xt4*logU3Q3*(
        pow4(mQ3)*(120*msq2*mU32 + 4241*pow4(mU3)) - 3459*mU32*pow6(mQ3) + 5*(
        24*msq2 - 571*mU32)*pow6(mU3) + mQ32*(-240*msq2*pow4(mU3) + 2039*pow6(
        mU3)) + 66*pow8(mQ3)))/(mU32*pow5(mQ32 - mU32)) + (16*lmQ3MR*pow2(Xt)*(
        (-60*logmqQ3)/mU32 - (1152*mQ32*dilog(1 - mQ32/mU32))/pow2(mQ32
        - mU32) - (2*(953*mQ32 - 377*mU32)*logU3Q3)/pow2(mQ32 - mU32) +
        (6*(36*mQ32 - 205*mU32))/(-(mQ32*mU32) + pow4(mU3)) - (18*lmQ3MR*logU3Q3*
        (-170*mQ32*mU32 + 37*pow4(mQ3) + 69*pow4(mU3)))/pow3(mQ32 -
        mU32) + (3*logU3Q3*(-1465*mU32*pow4(mQ3) + 4098*mQ32*pow4(mU3) +
        66*pow6(mQ3) - 1147*pow6(mU3)))/(mU32*pow3(-mQ32 + mU32)) - (6*pow2(
        logU3Q3)*(-383*mU32*pow4(mQ3) + 562*mQ32*pow4(mU3) + 101*pow6(
        mQ3) - 276*pow6(mU3)))/pow4(mQ32 - mU32) + (6*lmQ3MR*(-163*mU32*pow4(
        mQ3) + 470*mQ32*pow4(mU3) + 16*pow6(mQ3) - 131*pow6(mU3)))/(pow2(mQ32 -
        mU32)*pow4(mU3)) - (6*(6*pow4(mQ3)*(5*msq2*mU32 + 217*pow4(mU3)) - 321*
        mU32*pow6(mQ3) + 3*(10*msq2 + mU32)*pow6(mU3) - 4*mQ32*(15*msq2*pow4(
        mU3) + 61*pow6(mU3)) + 32*pow8(mQ3)))/(pow2(-(mQ3*mU32) + pow3(mQ3))*
        pow4(mU3))))/3. + (16*lmQ3MR*Xt4*pow2(logU3Q3)*(1325*pow4(mQ3)*
        pow4(mU3) - 746*mU32*pow6(mQ3) - 306*mQ32*pow6(mU3) + 80*pow8(mQ3) -
        345*pow8(mU3)))/pow6(mQ32 - mU32) + (8*lmQ3MR*logU3Q3*(3551*
        pow4(mQ3)*pow4(mU3) - 1775*mU32*pow6(mQ3) - 2679*mQ32*pow6(mU3) + 66*
        pow8(mQ3) + 877*pow8(mU3)))/(mU32*pow3(mQ32 - mU32)) + (16*lmQ3MR*Xt4*(
        5*(6*msq2*mU32 - 211*pow4(mU3))*pow6(mQ3) + pow4(mQ3)*(-30*msq2*pow4(
        mU3) + 3014*pow6(mU3)) - 417*mU32*pow8(mQ3) + 3*(10*msq2 + mU32)*pow8(
        mU3) - 3*mQ32*(10*msq2*pow6(mU3) + 523*pow8(mU3)) + 32*power10(mQ3)))/(
        mQ32*pow4(mU3)*pow4(mQ32 - mU32)) + (8*lmQ3MR*(3*(20*msq2*mU32 + 501*
        pow4(mU3))*pow6(mQ3) - 4*pow4(mQ3)*(15*msq2*pow4(mU3) + 259*pow6(mU3))
        - 706*mU32*pow8(mQ3) + 6*(10*msq2 + mU32)*pow8(mU3) + mQ32*(-60*msq2*
        pow6(mU3) + 177*pow8(mU3)) + 64*power10(mQ3)))/(pow2(-(mQ3*mU32) +
        pow3(mQ3))*pow4(mU3)))/3.;
      }

      case Limits::MU3_EQ_M3: {
         const double logmqQ3 = std::log(msq2/mQ32);
         const double logU3Q3 = std::log(mU32/mQ32);

         return (3*dlatas2Const - (32*lmQ3MR*(2701*mQ32 - 3853*mU32))/(9.*(mQ32 - mU32))
        + (8*lmQ3MR*(3857 + (216*mU32)/mQ32))/3. - 384*lmQ3MR*z3 + 160*lmQ3MR*
        logmqQ3 + (160*lmQ3MR*mU32*logmqQ3)/mQ32 - 4000*pow2(
        lmQ3MR) + (2352*mU32*pow2(lmQ3MR))/mQ32 + (16*(207*mQ32 - 271*mU32)*
        logU3Q3*pow2(lmQ3MR))/(mQ32 - mU32) + (lmQ3MR*(3321*mQ32 - 6137*
        mU32)*pow2(logU3Q3))/(mQ32 - mU32) + (160*lmQ3MR*(mQ32 + mU32)*
        Xt4*logmqQ3)/pow2(-(mQ3*mU32) + pow3(mQ3)) + 2208*pow3(lmQ3MR) +
        (32*lmQ3MR*(69*mQ32 - 73*mU32)*Xt4*dilog(1 - mU32/mQ32))/pow3(mQ32 -
        mU32) + (320*lmQ3MR*mU32*Xt4*logmqQ3*logU3Q3)/pow3(mQ32 -
        mU32) + (48*(69*mQ32 + 167*mU32)*Xt4*logU3Q3*pow2(lmQ3MR))/pow3(
        mQ32 - mU32) + (16*Xt4*pow2(lmQ3MR)*(179*mQ32*mU32 + 545*pow4(mQ3) -
        16*pow4(mU3)))/(pow2(mQ32 - mU32)*pow4(mQ3)) + (64*lmQ3MR*(-41*mQ32 +
        37*mU32)*dilog(1 - mU32/mQ32)*pow4(mU3))/pow3(-mQ32 + mU32) - (256*
        pow2(lmQ3MR)*pow4(mU3))/pow4(mQ3) - (64*lmQ3MR*mU3*(-5*mQ32 + mU32)*
        Xt5)/pow4(mQ32 - mU32) - (8*lmQ3MR*mU3*(313*mQ32 + 633*mU32)*Xt5*pow2(
        logU3Q3))/pow4(mQ32 - mU32) - (lmQ3MR*Xt4*pow2(logU3Q3)*(
        -480*mQ32*mU32 + 9*pow4(mQ3) + 5719*pow4(mU3)))/pow4(mQ32 - mU32) - (
        lmQ3MR*logU3Q3*(16*(377*mQ32 - 569*mU32)*pow2(mQ32 - mU32) + 3*
        logU3Q3*(-5131*mU32*pow4(mQ3) + 6267*mQ32*pow4(mU3) + 1113*pow6(
        mQ3) - 2313*pow6(mU3))))/(3.*pow3(mQ32 - mU32)) + (64*lmQ3MR*mU3*Xt*(2*
        ((-183 + 16*lmQ3MR)*mQ32 + 8*(5 - 2*lmQ3MR)*mU32)*pow2(mQ32 - mU32) +
        8*dilog(1 - mU32/mQ32)*(48*mU32*pow4(mQ3) - 47*mQ32*pow4(mU3)) + pow2(
        logU3Q3)*(-182*mU32*pow4(mQ3) + 77*mQ32*pow4(mU3) + 123*pow6(
        mQ3)) + logU3Q3*(10*logmqQ3*pow2(-(mQ3*mU32) + pow3(mQ3))
        + 12*(41 - 23*lmQ3MR)*mU32*pow4(mQ3) + (-199 + 138*lmQ3MR)*mQ32*pow4(
        mU3) + (-221 + 138*lmQ3MR)*pow6(mQ3) - 64*pow6(mU3))))/(mQ32*pow3(mQ32
        - mU32)) - (16*lmQ3MR*mU3*Xt5*logU3Q3*(-714*mU32*pow4(mQ3) +
        235*mQ32*pow4(mU3) + 383*pow6(mQ3) + 80*pow6(mU3)))/(mQ32*pow5(mQ32 -
        mU32)) - (lmQ3MR*logU3Q3*(-24902*mU32*pow4(mQ3) + 13707*mQ32*
        pow4(mU3) + 10779*pow6(mQ3) + 352*pow6(mU3)))/pow2(-(mQ3*mU32) + pow3(
        mQ3)) + (16*lmQ3MR*mU3*Xt5*(4*(201*mQ32 + 16*mU32)*pow2(mQ32 - mU32) +
        logU3Q3*(-1186*mU32*pow4(mQ3) + 41*mQ32*pow4(mU3) + 1081*pow6(
        mQ3) + 80*pow6(mU3))))/(mQ32*pow5(mQ32 - mU32)) + (8*lmQ3MR*mU3*Xt5*
        logU3Q3*(4*(185*mQ32 + 249*mU32)*pow2(mQ32 - mU32) + logU3Q3*(407*mU32*
        pow4(mQ3) - 1817*mQ32*pow4(mU3) + 937*pow6(mQ3) + 505*
        pow6(mU3))))/pow6(mQ32 - mU32) + (lmQ3MR*Xt4*logU3Q3*(16*pow2(
        mQ32 - mU32)*(-744*mQ32*mU32 + 377*pow4(mQ3) - 1553*pow4(mU3)) + 3*
        logU3Q3*(624*pow4(mQ3)*pow4(mU3) - 9554*mU32*pow6(mQ3) + 16658*mQ32*
        pow6(mU3) + 1113*pow8(mQ3) - 8713*pow8(mU3))))/(3.*pow6(mQ32 - mU32)) -
        (2*lmQ3MR*Xt4*logU3Q3*(3*pow4(mQ3)*(320*msq2*mU32 - 6251*pow4(
        mU3)) + (-480*msq2 + 769*mU32)*pow6(mQ3) + mQ32*(-480*msq2*pow4(mU3) +
        12379*pow6(mU3)) + 5205*pow8(mQ3) + 336*pow8(mU3)))/(mQ32*pow5(mQ32 -
        mU32)) + (64*lmQ3MR*mU3*pow3(Xt)*(4*(31*mQ32 - 33*mU32)*dilog(1 - mU32/
        mQ32)*pow2(-(mQ3*mU32) + pow3(mQ3)) - 10*logmqQ3*(2*(mQ32 -
        mU32) + (mQ32 + mU32)*logU3Q3)*pow2(-(mQ3*mU32) + pow3(mQ3)) +
        2248*pow4(mQ3)*pow4(mU3) - 828*lmQ3MR*pow4(mQ3)*pow4(mU3) - 2230*logU3Q3*
        pow4(mQ3)*pow4(mU3) + 266*lmQ3MR*logU3Q3*pow4(mQ3)*
        pow4(mU3) + 515*pow2(logU3Q3)*pow4(mQ3)*pow4(mU3) - 2408*mU32*
        pow6(mQ3) + 892*lmQ3MR*mU32*pow6(mQ3) + 1274*mU32*logU3Q3*pow6(
        mQ3) + 74*lmQ3MR*mU32*logU3Q3*pow6(mQ3) + 191*mU32*pow2(logU3Q3)*pow6(mQ3)
        - 568*mQ32*pow6(mU3) + 212*lmQ3MR*mQ32*pow6(mU3)
        + 998*mQ32*logU3Q3*pow6(mU3) - 202*lmQ3MR*mQ32*logU3Q3*
        pow6(mU3) - 475*mQ32*pow2(logU3Q3)*pow6(mU3) + 824*pow8(mQ3) -
        308*lmQ3MR*pow8(mQ3) - 106*logU3Q3*pow8(mQ3) - 138*lmQ3MR*logU3Q3*pow8(mQ3)
        - 263*pow2(logU3Q3)*pow8(mQ3) - 96*pow8(
        mU3) + 32*lmQ3MR*pow8(mU3) + 64*logU3Q3*pow8(mU3)))/(mQ32*pow5(
        mQ32 - mU32)) + (16*lmQ3MR*Xt4*((-30*msq2*mU32 + 3014*pow4(mU3))*pow6(
        mQ3) - 5*pow4(mQ3)*(6*msq2*pow4(mU3) + 211*pow6(mU3)) + 3*(10*msq2 -
        523*mU32)*pow8(mQ3) + mQ32*(30*msq2*pow6(mU3) - 417*pow8(mU3)) + 3*
        power10(mQ3) + 32*power10(mU3)))/(mU32*pow4(-(mQ3*mU32) + pow3(mQ3))) +
        (8*lmQ3MR*(-4*(15*msq2*mU32 + 259*pow4(mU3))*pow6(mQ3) + pow4(mQ3)*(-
        60*msq2*pow4(mU3) + 1503*pow6(mU3)) + 3*(20*msq2 + 59*mU32)*pow8(mQ3) +
        mQ32*(60*msq2*pow6(mU3) - 706*pow8(mU3)) + 6*power10(mQ3) + 64*power10(
        mU3)))/(mU32*pow2(mQ32 - mU32)*pow4(mQ3)) + (lmQ3MR*logU3Q3*(
        28393*pow4(mU3)*pow6(mQ3) - 4227*pow4(mQ3)*pow6(mU3) - 28089*mU32*pow8(
        mQ3) - 5904*mQ32*pow8(mU3) + 9635*power10(mQ3) + 512*power10(mU3)))/(
        pow3(mQ32 - mU32)*pow4(mQ3)) - (16*lmQ3MR*pow2(Xt)*(18*pow12(mQ3) +
        192*pow12(mU3) - 96*lmQ3MR*pow12(mU3) - 192*logU3Q3*pow12(mU3) +
        1152*dilog(1 - mU32/mQ32)*pow2(mQ32 - mU32)*pow4(mQ3)*pow4(mU3) + 60*
        mQ32*mU32*logmqQ3*pow4(mQ32 - mU32) + 1080*msq2*pow4(mU3)*pow6(
        mQ3) - 720*msq2*pow4(mQ3)*pow6(mU3) - 14676*pow6(mQ3)*pow6(mU3) + 7404*
        lmQ3MR*pow6(mQ3)*pow6(mU3) + 2325*logU3Q3*pow6(mQ3)*pow6(mU3) +
        3726*lmQ3MR*logU3Q3*pow6(mQ3)*pow6(mU3) + 5154*pow2(logU3Q3)*pow6(mQ3)*
        pow6(mU3) - 720*msq2*mU32*pow8(mQ3) + 6852*pow4(mU3)*
        pow8(mQ3) - 4392*lmQ3MR*pow4(mU3)*pow8(mQ3) + 3777*logU3Q3*pow4(
        mU3)*pow8(mQ3) - 4302*lmQ3MR*logU3Q3*pow4(mU3)*pow8(mQ3) - 5232*
        pow2(logU3Q3)*pow4(mU3)*pow8(mQ3) + 180*mQ32*msq2*pow8(mU3) +
        9978*pow4(mQ3)*pow8(mU3) - 4872*lmQ3MR*pow4(mQ3)*pow8(mU3) - 6817*logU3Q3*
        pow4(mQ3)*pow8(mU3) - 666*lmQ3MR*logU3Q3*pow4(mQ3)*
        pow8(mU3) - 726*pow2(logU3Q3)*pow4(mQ3)*pow8(mU3) + 180*msq2*
        power10(mQ3) - 270*mU32*power10(mQ3) + 786*lmQ3MR*mU32*power10(mQ3) -
        1175*mU32*logU3Q3*power10(mQ3) + 1242*lmQ3MR*mU32*logU3Q3
        *power10(mQ3) + 828*mU32*pow2(logU3Q3)*power10(mQ3) - 2094*mQ32*
        power10(mU3) + 1170*lmQ3MR*mQ32*power10(mU3) + 2082*mQ32*logU3Q3
        *power10(mU3)))/(3.*mU32*pow4(-(mQ3*mU32) + pow3(mQ3))) + (2*lmQ3MR*
        Xt4*(16*pow2(-(mQ3*mU32) + pow3(mQ3))*(-1343*mQ32*mU32 + 437*pow4(mQ3)
        - 54*pow4(mU3)) + 3*logU3Q3*(15523*pow4(mU3)*pow6(mQ3) - 977*
        pow4(mQ3)*pow6(mU3) - 14211*mU32*pow8(mQ3) - 2952*mQ32*pow8(mU3) +
        2425*power10(mQ3) + 256*power10(mU3))))/(3.*pow4(mQ3)*pow5(mQ32 - mU32)
        ))/3.;
      }

      case Limits::MQ3_EQ_MU3_EQ_M3: {
         const double logmqQ3 = std::log(msq2/mQ32);

         return dlatas2Const + (160*lmQ3MR*logmqQ3*(Xt4 - 6*mQ32*pow2(Xt) - 12*Xt*
        pow3(mQ3) + 2*mQ3*pow3(Xt) + 6*pow4(mQ3)))/(9.*pow4(mQ3)) + (8*lmQ3MR*(
        180*msq2*Xt4 + 176*mQ3*Xt5 - mQ32*((98 + 27*lmQ3MR)*Xt4 + 1080*msq2*
        pow2(Xt)) + 8*(-416 + 111*lmQ3MR)*pow3(mQ3)*pow3(Xt) + 24*(45*msq2 + 2*
        (8 + 57*lmQ3MR)*pow2(Xt))*pow4(mQ3) + 16*(100 - 477*lmQ3MR)*Xt*pow5(
        mQ3) + (-851 - 990*lmQ3MR - 432*z3 + 2484*pow2(lmQ3MR))*pow6(mQ3)))/(
        27.*pow6(mQ3));
      }

      default:
         break;
   }

   throw std::runtime_error("Mass limit not included!");
}

/**
 * Returns the shift needed to convert the 3L threshold correction of lambda to the MSbar scheme
 * @param xtOrder an integer key to omit the Xt contributions starting at xtOrder + 1
 * @param omitLogs an integer key to omit all log mu terms
 * @param omitXtLogs an integer key to omit all Xt^4*Log[mu] and Xt^5*Log[mu] terms
 */
double ThresholdCalculator::getDRbarPrimeToMSbarShift(int xtOrder, int omitLogs, int omitXtLogs) const
{
   const auto limit = static_cast<Limits>(p.massLimit3LThreshold);

   const double xtTerms =
      xtOrder <= 3
      ?
      getDRbarPrimeToMSbarXtTerms(limit, 4, omitXtLogs) +
      getDRbarPrimeToMSbarXtTerms(limit, 5, omitXtLogs) +
      getDRbarPrimeToMSbarXtTerms(limit, 6, omitXtLogs)
      :
      getDRbarPrimeToMSbarXtTerms(limit, 5, omitXtLogs) +
      getDRbarPrimeToMSbarXtTerms(limit, 6, omitXtLogs);

   const double g3as = getDeltaG3Alphas(omitLogs);
   const double ytas = getDeltaYtAlphas(limit, omitLogs);
   const double ytas2 = getDeltaYtAlphas2(limit, omitLogs);
   const double lambdaat = getDeltaLambdaAlphat(limit, omitLogs);
   const double lambdaatas = getDeltaLambdaAlphatAlphas(limit, omitLogs);

   return -(-2.*(lambdaat*(3*pow2(ytas) + 2*ytas2) + (lambdaatas - 4*ytas*lambdaat)
      *(g3as + 2*ytas)) - xtTerms);
}

/**
 * Returns the DRbarPrime to MSbar shift of delta lambda 3L at a given xtOrder
 * @param limit an integer key for a mass limit
 * @param xtOrder an integer key to omit the Xt contributions starting at xtOrder + 1
 * @param omitLogs an integer key to omit all log mu terms
 */
double ThresholdCalculator::getDRbarPrimeToMSbarXtTerms(Limits limit, int xtOrder, int omitLogs) const
{
   using himalaya::dilog;

   const double Mst12 = pow2(p.MSt(0));
   const double m32 = pow2(p.MG);
   const double MR2 = pow2(p.scale);
   const double mU32 = p.mu2(2,2);
   const double mQ32 = p.mq2(2,2);
   const double m3 = std::sqrt(m32);
   const double mU3 = std::sqrt(mU32);
   const double mQ3 = std::sqrt(mQ32);
   const double msq = std::sqrt(msq2);
   const double Xt = p.Au(2,2) - p.mu * p.vd / p.vu;
   const double Xt4 = xtOrder != 4 ? 0 : pow4(Xt);
   const double Xt5 = xtOrder != 5 ? 0 : pow5(Xt);
   const double Xt6 = xtOrder != 6 ? 0 : pow6(Xt);
   const double lmQ3MR = omitLogs * std::log(Mst12 / MR2) + std::log(mQ32 / Mst12);

   switch (limit) {
      case Limits::GENERAL: {
         const double logQ3U3 = std::log(mQ32/mU32);
         const double logQ3m3 = std::log(mQ32/m32);
         const double logQ3mq = std::log(mQ32/msq2);

         return
         (-4*Xt4*(640*mQ32*(mQ32 - mU32)*mU32*logQ3U3*pow2(m32 - mQ32)*
         pow2(m32 - mU32)*pow2(m3*(m32 - mQ32)*mU32*logQ3U3 + (mQ32 -
         mU32)*logQ3m3*pow3(m3)) - 512*m3*mQ32*mU32*(-8*mQ32 + 4*lmQ3MR* mQ32 +
         8*mU32 - 4*lmQ3MR*mU32 + (4*m32 - 2*(mQ32 + mU32))*dilog(1 - m32/mQ32)
         + 2*(-2*m32 + mQ32 + mU32)*dilog(1 - m32/mU32) + 2*mQ32* logQ3U3 -
         2*lmQ3MR*mQ32*logQ3U3 + 6*mU32*logQ3U3 - 2* lmQ3MR*mU32*logQ3U3 +
         4*m32*logQ3m3*logQ3U3 - 2* m32*pow2(logQ3U3) + mQ32*pow2(logQ3U3) +
         mU32*pow2(logQ3U3))*(m3*(m32 - mQ32)*mU32*logQ3U3 + (mQ32 -
         mU32)*logQ3m3*pow3(m3))*pow3(m32 - mQ32)*pow3(m32 - mU32) - 4*pow2(m32
         - mQ32)*pow2(m32 - mU32)*(6*pow2(m32 - mQ32)*pow2(m32 - mU32) - 12*
         lmQ3MR*pow2(m32 - mQ32)*pow2(m32 - mU32) + 12*logQ3m3*pow2(m32 -
         mQ32)*pow2(m32 - mU32) - (12*lmQ3MR - 10*logQ3mq - logQ3U3)*pow2(m32 -
         mQ32)*pow2(m32 - mU32) + 16*(mQ32*(-m32 + mQ32)*(m32 - mU32)*mU32 +
         (m32 - mQ32)*(m32 - mU32)*pow4(m3) + lmQ3MR*pow2(m32 -
         mU32)*(2*m32*mQ32 - pow4(mQ3)) + (lmQ3MR - logQ3U3)*pow2(m32 -
         mQ32)*(2*m32*mU32 - pow4(mU3)) - (lmQ3MR - logQ3m3)*pow4(m3)*(-2*
         m32*(mQ32 + mU32) + 2*pow4(m3) + pow4(mQ3) + pow4(mU3))))*(mQ32*(mQ32
         - mU32)*mU32*(-2*m32 + mQ32 + mU32)*dilog(1 - m32/mQ32) +
         m32*mU32*pow4( mQ3) - lmQ3MR*m32*mU32*pow4(mQ3) +
         m32*mU32*logQ3m3*pow4(mQ3) + 2*m32*mU32*logQ3U3*pow4(mQ3) -
         2*lmQ3MR*m32*mU32*logQ3U3* pow4(mQ3) +
         m32*mU32*pow2(logQ3U3)*pow4(mQ3) + m32*mQ32*pow4( mU3) -
         lmQ3MR*m32*mQ32*pow4(mU3) + m32*mQ32*logQ3m3*pow4(mU3) -
         2*m32*mQ32*logQ3U3*pow4(mU3) + 2*lmQ3MR*m32*mQ32*logQ3U3* pow4(mU3) -
         m32*mQ32*pow2(logQ3U3)*pow4(mU3) + 12*pow4(mQ3)* pow4(mU3) -
         8*lmQ3MR*pow4(mQ3)*pow4(mU3) + 4*logQ3U3*pow4(mQ3)* pow4(mU3) -
         2*pow2(logQ3U3)*pow4(mQ3)*pow4(mU3) + mQ32*mU32* dilog(1 -
         m32/mU32)*(2*m32*(mQ32 - mU32) - pow4(mQ3) + pow4(mU3)) -
         m32*pow6(mQ3) + lmQ3MR*m32*pow6(mQ3) - 6*mU32*pow6(mQ3) + 4*lmQ3MR*
         mU32*pow6(mQ3) - m32*logQ3m3*pow6(mQ3) + 3*mU32*logQ3U3* pow6(mQ3) -
         2*lmQ3MR*mU32*logQ3U3*pow6(mQ3) - m32*pow6(mU3) + lmQ3MR*m32*pow6(mU3)
         - 6*mQ32*pow6(mU3) + 4*lmQ3MR*mQ32*pow6(mU3) - m32*logQ3m3*pow6(mU3) -
         7*mQ32*logQ3U3*pow6(mU3) + 2* lmQ3MR*mQ32*logQ3U3*pow6(mU3) -
         2*mQ32*pow2(logQ3U3)* pow6(mU3)) + mQ32*mU32*(2*(mQ32 - mU32) - (mQ32
         + mU32)*logQ3U3) *(12*(mQ32 - mU32)*pow2(-((m32 - mQ32)*((m32 -
         mQ32)*(2*m32 - mU32)* mU32*logQ3U3 + (m32 - mU32)*((1 +
         2*lmQ3MR)*mQ32*mU32 - 2* lmQ3MR*m32*(mQ32 + mU32) + (-1 +
         2*lmQ3MR)*pow4(m3)))) + logQ3m3* pow4(m3)*(-2*m32*(mQ32 + mU32) +
         2*pow4(m3) + pow4(mQ3) + pow4(mU3))) - 15*(m32 - mQ32)*(m32 -
         mU32)*(mQ32 - mU32)*pow2(lmQ3MR - logQ3mq )*(-2*(mQ32 - msq2)*pow3(m32
         - mU32)*(-3*m32*(mQ32 - msq2) + mQ32*(mQ32 + msq2) - 2*pow4(m3)) -
         2*(msq2 - mU32)*pow3(m32 - mQ32)*(-3*m32*(msq2 - mU32) - mU32*(msq2 +
         mU32) + 2*pow4(m3))) + 60*(m32 - mU32)*(mQ32 - mU32)*(msq2 -
         mU32)*dilog(1 - msq2/mU32)*(-3*m32*(msq2 - mU32) - mU32*( msq2 + mU32)
         + 2*pow4(m3))*pow4(m32 - mQ32) - 60*(m32 - mQ32)*(mQ32 - msq2)*(mQ32 -
         mU32)*dilog(1 - msq2/mQ32)*(3*m32*(mQ32 - msq2) - mQ32*( mQ32 + msq2)
         + 2*pow4(m3))*pow4(m32 - mU32) - 20*(mQ32 - mU32)*(lmQ3MR -
         logQ3mq)*pow2(m32 - mQ32)*pow2(m32 - mU32)*(3*mQ32*msq2*pow4( mU3) +
         pow4(mQ3)*(3*msq2*mU32 + 10*pow4(mU3)) + pow4(m3)*(-15*msq2*mU32 +
         mQ32*(-15*msq2 + 44*mU32) + 11*pow4(mQ3) + 11*pow4(mU3)) + 3*m32*((3*
         msq2 - 7*mU32)*pow4(mQ3) + 3*msq2*pow4(mU3) - mQ32*(4*msq2*mU32 + 7*
         pow4(mU3))) + (-23*mQ32 + 18*msq2 - 23*mU32)*pow6(m3) + 12*pow8(m3)) +
         (mQ32 - mU32)*pow2(m32 - mQ32)*pow2(m32 - mU32)*(-(mQ32*mU32*(6*mU32*(
         10*msq2 + mU32) + mQ32*(60*msq2 + 169*mU32) + 6*pow4(mQ3))) +
         pow4(m3)* (4*mQ32*(75*msq2 - 28*mU32) + 300*msq2*mU32 + 93*pow4(mQ3) +
         93*pow4( mU3)) - 6*(47*mQ32 + 60*msq2 + 47*mU32)*pow6(m3) -
         2*m32*(2*(45*msq2 - 52*mU32)*pow4(mQ3) + 9*(10*msq2 + mU32)*pow4(mU3)
         - 8*mQ32*(15*msq2* mU32 + 13*pow4(mU3)) + 9*pow6(mQ3)) + 291*pow8(m3))
         + 60*(m32 - mQ32)*( m32 - msq2)*(m32 - mU32)*(mQ32 - mU32)*dilog(1 -
         m32/msq2)*(-(mQ32* msq2*mU32*(pow4(mQ3) + pow4(mU3))) + (-8*msq2*mU32
         + mQ32*(-8*msq2 + 30*mU32) + 3*pow4(mQ3) + 3*pow4(mU3))*pow6(m3) -
         pow4(m3)*((-9*msq2 + 15*mU32)*pow4(mQ3) - 9*msq2*pow4(mU3) +
         3*mQ32*(2*msq2*mU32 + 5*pow4( mU3)) + pow6(mQ3) + pow6(mU3)) +
         m32*(3*msq2*mU32*pow4(mQ3) + (-3*msq2 + 5*mU32)*pow6(mQ3) -
         3*msq2*pow6(mU3) + mQ32*(3*msq2*pow4(mU3) + 5* pow6(mU3))) + (-8*mQ32
         + 6*msq2 - 8*mU32)*pow8(m3) + 2*power10(m3)) + 2*(m32 - mQ32)*(m32 -
         mU32)*dilog(1 - mQ32/mU32)*pow2(mQ32 - mU32)*((- 76*mQ32*mU32 +
         34*pow4(mQ3) + 34*pow4(mU3))*pow6(m3) + 7*pow4(mU3)* pow6(mQ3) +
         pow4(m3)*(65*mU32*pow4(mQ3) + 65*mQ32*pow4(mU3) - 29*pow6( mQ3) -
         29*pow6(mU3)) + 7*pow4(mQ3)*pow6(mU3) - 14*(mQ32 + mU32)*pow8( m3) +
         3*mU32*pow8(mQ3) + 3*mQ32*pow8(mU3) + m32*(-34*pow4(mQ3)*pow4( mU3) -
         26*mU32*pow6(mQ3) - 26*mQ32*pow6(mU3) + 9*pow8(mQ3) + 9*pow8( mU3)) +
         12*power10(m3)) - pow2(lmQ3MR - logQ3U3)*pow4(m32 - mQ32)*((-mQ32 +
         mU32)*pow4(mU3)*(4*mQ32*mU32 + 3*pow4(mQ3) + 30*pow4( msq) +
         21*pow4(mU3)) + (60*msq2*mU32 + mQ32*(-60*msq2 + 68*mU32) - 2*
         pow4(mQ3) + 702*pow4(mU3))*pow6(m3) + 2*m32*mU32*(18*mU32*pow4(mQ3) +
         6*mU32*(-15*msq2*mU32 + 5*pow4(msq) + 4*pow4(mU3)) +
         mQ32*(90*msq2*mU32 - 30*pow4(msq) + 25*pow4(mU3)) - 3*pow6(mQ3)) +
         pow4(m3)*(-33*mU32* pow4(mQ3) - 90*mU32*pow4(msq) +
         3*mQ32*(-40*msq2*mU32 + 30*pow4(msq) - 47*pow4(mU3)) +
         120*msq2*pow4(mU3) + 9*pow6(mQ3) - 347*pow6(mU3)) + ( 44*mQ32 -
         556*mU32)*pow8(m3) + 128*power10(m3)) - 2*(mQ32 - mU32)*pow2( lmQ3MR -
         logQ3m3)*(878*(mQ32 + mU32)*pow14(m3) - 262*pow16(m3) +
         pow4(m3)*pow4(mQ3)*pow4(mU3)*(452*mQ32*mU32 + pow4(mQ3) + pow4(mU3)) -
         pow12(m3)*(2972*mQ32*mU32 + 885*pow4(mQ3) + 885*pow4(mU3)) - 144*m32*(
         mQ32 + mU32)*pow6(mQ3)*pow6(mU3) + 36*pow8(mQ3)*pow8(mU3) - pow8(m3)*(
         2404*pow4(mQ3)*pow4(mU3) + 1004*mU32*pow6(mQ3) + 1004*mQ32*pow6(mU3) +
         49*pow8(mQ3) + 49*pow8(mU3)) + 2*pow6(m3)*(184*pow4(mU3)*pow6(mQ3) +
         184*pow4(mQ3)*pow6(mU3) + 79*mU32*pow8(mQ3) + 79*mQ32*pow8(mU3)) + 4*(
         733*mU32*pow4(mQ3) + 733*mQ32*pow4(mU3) + 80*pow6(mQ3) +
         80*pow6(mU3))* power10(m3)) + 2*lmQ3MR*(m32 - mQ32)*(m32 -
         mU32)*(-10*(mQ32 - mU32)*( lmQ3MR - logQ3mq)*pow3(m32 - mU32)*(2*(mQ32
         + 3*msq2)*pow4(m3) - 3*mQ32*pow4(msq) - 3*m32*(-6*mQ32*msq2 +
         pow4(mQ3) + 3*pow4(msq)) + pow6(mQ3)) + (m32 - mU32)*((mQ32 -
         mU32)*mU32*pow4(mQ3)*(128*mQ32*mU32 + 3*pow4(mQ3) - 3*(10*msq2*mU32 +
         pow4(mU3))) + mQ32*pow4(m3)*(-5*(6* msq2 + 5*mU32)*pow4(mQ3) +
         6*mQ32*(35*msq2*mU32 - 157*pow4(mU3)) - 180* msq2*pow4(mU3) +
         106*pow6(mQ3) - 163*pow6(mU3)) + pow6(m3)*((-90*msq2 +
         985*mU32)*pow4(mQ3) + mQ32*(90*msq2*mU32 + 622*pow4(mU3)) - 18*pow6(
         mQ3) - 53*pow6(mU3)) + (-764*mQ32*mU32 - 365*pow4(mQ3) +
         105*pow4(mU3)) *pow8(m3) + m32*mQ32*(6*pow4(mQ3)*(10*msq2*mU32 +
         41*pow4(mU3)) - 279* mU32*pow6(mQ3) + 9*(10*msq2 + mU32)*pow6(mU3) +
         mQ32*(-150*msq2*pow4( mU3) + 271*pow6(mU3)) + 9*pow8(mQ3)) +
         4*(77*mQ32 - 13*mU32)*power10( m3)) - (mQ32 - mU32)*(lmQ3MR -
         logQ3U3)*(mU32*(4*mQ32*mU32 + 3* pow4(mQ3) - 3*pow4(mU3))*pow6(mQ3) +
         pow6(m3)*(-53*mU32*pow4(mQ3) + 57* mQ32*pow4(mU3) + 35*pow6(mQ3) +
         pow6(mU3)) + m32*pow4(mQ3)*(-35*mU32* pow4(mQ3) + mQ32*pow4(mU3) +
         9*pow6(mQ3) + 9*pow6(mU3)) - (28*mQ32*mU32 + 17*pow4(mQ3) +
         3*pow4(mU3))*pow8(m3) - pow4(m3)*(23*pow4(mQ3)*pow4( mU3) -
         75*mU32*pow6(mQ3) + 19*mQ32*pow6(mU3) + 29*pow8(mQ3)) + 2*(7* mQ32 +
         mU32)*power10(m3))) - 2*(lmQ3MR - logQ3m3)*(20*(m32 - mQ32)*(m32 -
         mU32)*(mQ32 - mU32)*(lmQ3MR - logQ3mq)*pow4(m3)*(7*
         mQ32*mU32*(pow4(mQ3) + pow4(mU3)) + 3*pow4(m3)*(14*mQ32*mU32 + pow4(
         mQ3) + pow4(mU3)) - 10*(mQ32 + mU32)*pow6(m3) - m32*(21*mU32*pow4(mQ3)
         + 21*mQ32*pow4(mU3) + pow6(mQ3) + pow6(mU3)) + 2*pow8(m3)) + (m32 -
         mQ32)*(m32 - mU32)*(mQ32 - mU32)*(864*pow12(m3) + 48*pow6(mQ3)*pow6(
         mU3) + 3*m32*mQ32*mU32*(5*(2*msq2 - 23*mU32)*pow4(mQ3) -
         115*mQ32*pow4( mU3) + 10*msq2*pow4(mU3) + pow6(mQ3) + pow6(mU3)) -
         9*pow6(m3)*((30* msq2 + 317*mU32)*pow4(mQ3) + 30*msq2*pow4(mU3) +
         mQ32*(-20*msq2*mU32 + 317*pow4(mU3)) + 39*pow6(mQ3) + 39*pow6(mU3)) +
         4*(60*msq2*mU32 + 2* mQ32*(30*msq2 + 511*mU32) + 365*pow4(mQ3) +
         365*pow4(mU3))*pow8(m3) + pow4(m3)*(pow4(mQ3)*(-90*msq2*mU32 +
         1894*pow4(mU3)) + (90*msq2 + 572* mU32)*pow6(mQ3) + 9*(10*msq2 +
         mU32)*pow6(mU3) + mQ32*(-90*msq2*pow4( mU3) + 572*pow6(mU3)) +
         9*pow8(mQ3)) - 2*(971*mQ32 + 90*msq2 + 971* mU32)*power10(m3)) -
         2*(m32 - mQ32)*(lmQ3MR - logQ3U3)*(-12*(17* mQ32 + 15*mU32)*pow14(m3)
         + 64*pow16(m3) + 2*pow12(m3)*(302*mQ32*mU32 + 117*pow4(mQ3) +
         61*pow4(mU3)) + 6*m32*pow4(mQ3)*(-(mQ32*mU32) + 4*pow4( mQ3) -
         3*pow4(mU3))*pow6(mU3) + mQ32*pow4(m3)*pow4(mU3)*(-53*mU32*pow4( mQ3)
         + 61*mQ32*pow4(mU3) + 43*pow6(mQ3) + 13*pow6(mU3)) + 6*(-mQ32 +
         mU32)*pow6(mQ3)*pow8(mU3) - mU32*pow6(m3)*(13*pow4(mQ3)*pow4(mU3) +
         239*mU32*pow6(mQ3) + 49*mQ32*pow6(mU3) + 76*pow8(mQ3) + 7*pow8(mU3)) +
         pow8(m3)*(503*pow4(mQ3)*pow4(mU3) + 401*mU32*pow6(mQ3) + 13*mQ32*pow6(
         mU3) + 13*pow8(mQ3) + 30*pow8(mU3)) - (779*mU32*pow4(mQ3) + 369*mQ32*
         pow4(mU3) + 101*pow6(mQ3) + 31*pow6(mU3))*power10(m3)) + 2*lmQ3MR*(m32
         - mU32)*(-12*(15*mQ32 + 17*mU32)*pow14(m3) + 64*pow16(m3) +
         2*pow12(m3) *(302*mQ32*mU32 + 61*pow4(mQ3) + 117*pow4(mU3)) -
         6*m32*(mQ32*mU32 + 3* pow4(mQ3) - 4*pow4(mU3))*pow4(mU3)*pow6(mQ3) +
         mU32*pow4(m3)*pow4(mQ3)* (61*mU32*pow4(mQ3) - 53*mQ32*pow4(mU3) +
         13*pow6(mQ3) + 43*pow6(mU3)) + 6*(mQ32 - mU32)*pow6(mU3)*pow8(mQ3) +
         pow8(m3)*(503*pow4(mQ3)*pow4(mU3) + 13*mU32*pow6(mQ3) +
         401*mQ32*pow6(mU3) + 30*pow8(mQ3) + 13*pow8(mU3)) -
         mQ32*pow6(m3)*(13*pow4(mQ3)*pow4(mU3) + 49*mU32*pow6(mQ3) + 239*mQ32*
         pow6(mU3) + 7*pow8(mQ3) + 76*pow8(mU3)) - (369*mU32*pow4(mQ3) + 779*
         mQ32*pow4(mU3) + 31*pow6(mQ3) + 101*pow6(mU3))*power10(m3))) + 2*(m32
         - mU32)*dilog(1 - m32/mQ32)*pow2(m32 - mQ32)*(128*pow12(m3) -
         mU32*pow4( mQ3)*(mU32*pow4(mQ3) + 19*mQ32*pow4(mU3) + 3*pow6(mQ3) -
         23*pow6(mU3)) + pow6(m3)*(-251*mU32*pow4(mQ3) - 1025*mQ32*pow4(mU3) +
         9*pow6(mQ3) - 13*pow6(mU3)) + (902*mQ32*mU32 + 98*pow4(mQ3) +
         280*pow4(mU3))*pow8(m3) + pow4(m3)*(417*pow4(mQ3)*pow4(mU3) -
         151*mU32*pow6(mQ3) + 391*mQ32* pow6(mU3) + 20*pow8(mQ3) -
         37*pow8(mU3)) - 2*(147*mQ32 + 173*mU32)* power10(m3) +
         m32*(41*pow4(mU3)*pow6(mQ3) - 167*pow4(mQ3)*pow6(mU3) +
         41*mU32*pow8(mQ3) - 34*mQ32*pow8(mU3) - 9*power10(mQ3))) + (m32 -
         mU32) *pow2(lmQ3MR)*(-4*(139*mQ32 + 85*mU32)*pow14(m3) + 128*pow16(m3)
         + 2* pow12(m3)*(-30*msq2*mU32 + mQ32*(30*msq2 + 856*mU32) +
         357*pow4(mQ3) + 131*pow4(mU3)) - (mQ32 -
         mU32)*mU32*pow4(mQ3)*(21*pow4(mQ3)*pow4(mU3) + 30*pow4(msq)*pow4(mU3)
         + 4*mU32*pow6(mQ3) + 3*pow8(mQ3)) + pow8(m3)*(-(
         pow4(mU3)*(180*msq2*mU32 + 270*pow4(msq) + 43*pow4(mU3))) +
         pow4(mQ3)*( -180*msq2*mU32 + 60*pow4(msq) + 2771*pow4(mU3)) +
         (-180*msq2 + 933* mU32)*pow6(mQ3) + 3*mQ32*(70*mU32*pow4(msq) +
         180*msq2*pow4(mU3) + 241* pow6(mU3)) + 96*pow8(mQ3)) + ((120*msq2 -
         2209*mU32)*pow4(mQ3) - mQ32*( 300*msq2*mU32 + 90*pow4(msq) +
         1903*pow4(mU3)) - 373*pow6(mQ3) + 5*(18* mU32*pow4(msq) +
         36*msq2*pow4(mU3) + pow6(mU3)))*power10(m3) - 2*pow6(
         m3)*(-3*(90*msq2*mU32 + 5*pow4(msq) - 251*pow4(mU3))*pow6(mQ3) + pow4(
         mQ3)*(105*mU32*pow4(msq) + 90*msq2*pow4(mU3) + 593*pow6(mU3)) - 53*
         mU32*pow8(mQ3) + 15*mQ32*(3*pow4(msq)*pow4(mU3) + 14*msq2*pow6(mU3) +
         2*pow8(mU3)) - 15*(9*pow4(msq)*pow6(mU3) + 2*msq2*pow8(mU3)) + 21*
         power10(mQ3)) + 2*pow4(m3)*(19*pow12(mQ3) - 5*pow6(mQ3)*(9*mU32*pow4(
         msq) + 54*msq2*pow4(mU3) - 61*pow6(mU3)) + 158*pow4(mU3)*pow8(mQ3) -
         45*pow4(msq)*pow8(mU3) + 3*pow4(mQ3)*(45*pow4(msq)*pow4(mU3) +
         70*msq2* pow6(mU3) + 27*pow8(mU3)) + mQ32*(-45*pow4(msq)*pow6(mU3) +
         60*msq2* pow8(mU3)) - 115*mU32*power10(mQ3)) + m32*(47*mU32*pow12(mQ3)
         - 9* pow14(mQ3) - 119*pow6(mU3)*pow8(mQ3) +
         6*pow6(mQ3)*(15*pow4(msq)*pow4( mU3) + 30*msq2*pow6(mU3) -
         14*pow8(mU3)) + 60*mQ32*pow4(msq)*pow8(mU3) -
         30*pow4(mQ3)*(5*pow4(msq)*pow6(mU3) + 6*msq2*pow8(mU3)) + 37*pow4(
         mU3)*power10(mQ3))) - 2*(m32 - mQ32)*dilog(1 - m32/mU32)*pow2(m32 -
         mU32)*(128*pow12(m3) - pow6(m3)*(1025*mU32*pow4(mQ3) + 251*mQ32*pow4(
         mU3) + 13*pow6(mQ3) - 9*pow6(mU3)) +
         mQ32*pow4(mU3)*(-19*mU32*pow4(mQ3) - mQ32*pow4(mU3) + 23*pow6(mQ3) -
         3*pow6(mU3)) + (902*mQ32*mU32 + 280* pow4(mQ3) +
         98*pow4(mU3))*pow8(m3) + pow4(m3)*(417*pow4(mQ3)*pow4(mU3) +
         391*mU32*pow6(mQ3) - 151*mQ32*pow6(mU3) - 37*pow8(mQ3) + 20*pow8(mU3)
         ) - 2*(173*mQ32 + 147*mU32)*power10(m3) +
         m32*(-167*pow4(mU3)*pow6(mQ3) + 41*pow4(mQ3)*pow6(mU3) -
         34*mU32*pow8(mQ3) + 41*mQ32*pow8(mU3) - 9* power10(mU3))) + 2*(m32 -
         mU32)*(lmQ3MR - logQ3U3)*pow2(m32 -
         mQ32)*(180*msq2*mU32*pow4(m3)*pow4(mQ3) - 210*mQ32*msq2*pow4(m3)*pow4(
         mU3) + 150*m32*msq2*pow4(mQ3)*pow4(mU3) + 942*pow4(m3)*pow4(mQ3)*pow4(
         mU3) - 90*mQ32*msq2*mU32*pow6(m3) - 622*mU32*pow4(mQ3)*pow6(m3) - 985*
         mQ32*pow4(mU3)*pow6(m3) + 90*msq2*pow4(mU3)*pow6(m3) - 90*m32*msq2*
         mU32*pow6(mQ3) + 163*mU32*pow4(m3)*pow6(mQ3) - 271*m32*pow4(mU3)*pow6(
         mQ3) - 30*msq2*pow4(mU3)*pow6(mQ3) + 53*pow6(m3)*pow6(mQ3) - 60*m32*
         mQ32*msq2*pow6(mU3) + 25*mQ32*pow4(m3)*pow6(mU3) + 30*msq2*pow4(m3)*
         pow6(mU3) - 246*m32*pow4(mQ3)*pow6(mU3) + 30*msq2*pow4(mQ3)*pow6(mU3)
         + 18*pow6(m3)*pow6(mU3) + 131*pow6(mQ3)*pow6(mU3) - 10*(mQ32 - mU32)*(
         lmQ3MR - logQ3mq)*pow2(m32 - mQ32)*(2*(3*msq2 + mU32)*pow4(m3) -
         3*mU32*pow4(msq) - 3*m32*(-6*msq2*mU32 + 3*pow4(msq) + pow4(mU3)) +
         pow6(mU3)) + 764*mQ32*mU32*pow8(m3) - 105*pow4(mQ3)*pow8(m3) + 365*
         pow4(mU3)*pow8(m3) - 9*m32*mU32*pow8(mQ3) - 3*pow4(mU3)*pow8(mQ3) +
         279*m32*mQ32*pow8(mU3) - 106*pow4(m3)*pow8(mU3) - 125*pow4(mQ3)*pow8(
         mU3) + 52*mQ32*power10(m3) - 308*mU32*power10(m3) - 9*m32*power10(mU3)
         - 3*mQ32*power10(mU3)))))/(3.*mQ32*mU32*pow4(m32 - mQ32)*pow4(m32 -
         mU32)*pow4(mQ32 - mU32))+ (-32*Xt5*(32*pow2(m32 - mQ32)*pow2(m32 -
         mU32)*(m3*(m32 - mQ32)*mU32* logQ3U3 + (mQ32 -
         mU32)*logQ3m3*pow3(m3))*(mQ32*(mQ32 - mU32)* mU32*(-2*m32 + mQ32 +
         mU32)*dilog(1 - m32/mQ32) + m32*mU32*pow4(mQ3) -
         lmQ3MR*m32*mU32*pow4(mQ3) + m32*mU32*logQ3m3*pow4(mQ3) + 2*m32*
         mU32*logQ3U3*pow4(mQ3) - 2*lmQ3MR*m32*mU32*logQ3U3*pow4( mQ3) +
         m32*mU32*pow2(logQ3U3)*pow4(mQ3) + m32*mQ32*pow4(mU3) -
         lmQ3MR*m32*mQ32*pow4(mU3) + m32*mQ32*logQ3m3*pow4(mU3) - 2*m32*
         mQ32*logQ3U3*pow4(mU3) + 2*lmQ3MR*m32*mQ32*logQ3U3*pow4( mU3) -
         m32*mQ32*pow2(logQ3U3)*pow4(mU3) + 12*pow4(mQ3)*pow4(mU3) -
         8*lmQ3MR*pow4(mQ3)*pow4(mU3) + 4*logQ3U3*pow4(mQ3)*pow4(mU3) -
         2*pow2(logQ3U3)*pow4(mQ3)*pow4(mU3) + mQ32*mU32*dilog(1 - m32/
         mU32)*(2*m32*(mQ32 - mU32) - pow4(mQ3) + pow4(mU3)) - m32*pow6(mQ3) +
         lmQ3MR*m32*pow6(mQ3) - 6*mU32*pow6(mQ3) + 4*lmQ3MR*mU32*pow6(mQ3) -
         m32*logQ3m3*pow6(mQ3) + 3*mU32*logQ3U3*pow6(mQ3) - 2*
         lmQ3MR*mU32*logQ3U3*pow6(mQ3) - m32*pow6(mU3) + lmQ3MR*m32*pow6( mU3)
         - 6*mQ32*pow6(mU3) + 4*lmQ3MR*mQ32*pow6(mU3) - m32*logQ3m3* pow6(mU3)
         - 7*mQ32*logQ3U3*pow6(mU3) + 2*lmQ3MR*mQ32*logQ3U3*pow6(mU3) -
         2*mQ32*pow2(logQ3U3)*pow6(mU3)) - mQ32*mU32*( 2*(mQ32 - mU32) - (mQ32
         + mU32)*logQ3U3)*(-60*m3*(m32 - mQ32)*( m32 - msq2)*(m32 -
         mU32)*(mQ32*(msq2 - 2*mU32) + msq2*mU32 + m32*(mQ32 - 2*msq2 +
         mU32))*dilog(1 - m32/msq2)*pow2(mQ32 - mU32) - 60*m3*msq2*( lmQ3MR -
         logQ3mq)*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2(mQ32 - mU32) -
         30*m3*(m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)*(pow2(mQ32 -
         msq2)*pow2(m32 - mU32) - pow2(m32 - mQ32)*pow2(msq2 - mU32))*pow2(
         lmQ3MR - logQ3mq) + 2*pow2(m32 - mQ32)*pow2(m32 - mU32)*pow2( mQ32 -
         mU32)*(-3*m3*(mQ32 + 10*msq2 + mU32) + 14*pow3(m3)) + 60*m3*(m32 -
         mU32)*(mQ32 - mU32)*dilog(1 - msq2/mU32)*pow2(msq2 - mU32)*pow3(m32 -
         mQ32) - 60*m3*(m32 - mQ32)*(mQ32 - mU32)*dilog(1 -
         msq2/mQ32)*pow2(mQ32 - msq2)*pow3(m32 - mU32) + 2*m3*(m32 - mQ32)*(m32
         - mU32)*dilog(1 - mQ32/mU32)*pow3(mQ32 - mU32)*(-(mQ32*mU32) -
         5*m32*(mQ32 + mU32) + 5* pow4(m3) + 3*pow4(mQ3) + 3*pow4(mU3)) -
         2*m3*(m32 - mQ32)*(mQ32 - mU32) *dilog(1 - m32/mU32)*pow3(m32 -
         mU32)*(-7*mQ32*mU32 + m32*(-37*mQ32 + mU32) + 18*pow4(m3) +
         22*pow4(mQ3) + 3*pow4(mU3)) + 2*m3*(m32 - mU32)*( mQ32 - mU32)*dilog(1
         - m32/mQ32)*pow3(m32 - mQ32)*(m32*(mQ32 - 37*mU32) - 7*mQ32*mU32 +
         18*pow4(m3) + 3*pow4(mQ3) + 22*pow4(mU3)) + 12*m3*(mQ32 -
         mU32)*(m32*(mQ32 - mU32)*logQ3m3 + (m32 - mQ32)*mU32*logQ3U3)*(-((m32
         - mQ32)*((m32 - mQ32)*(2*m32 - mU32)*mU32*logQ3U3 + (m32 - mU32)*((1 +
         2*lmQ3MR)*mQ32*mU32 - 2*lmQ3MR*m32*(mQ32 + mU32) + (-1 +
         2*lmQ3MR)*pow4(m3)))) + logQ3m3*pow4(m3)*(-2*m32*(mQ32 + mU32) +
         2*pow4(m3) + pow4(mQ3) + pow4(mU3))) - 2*m3*(m32 - mU32)*(mQ32 -
         mU32)*(lmQ3MR - logQ3U3)*pow2(m32 - mQ32)*(-53*m32*mQ32*mU32 -
         30*m32*msq2*mU32 + 30*mQ32*msq2*mU32 - 16*mQ32*pow4(m3) +
         55*mU32*pow4( m3) + 3*mU32*pow4(mQ3) - 61*m32*pow4(mU3) +
         53*mQ32*pow4(mU3) + 5*(m32 - mQ32)*(lmQ3MR - logQ3mq)*(-(m32*mU32) -
         12*msq2*mU32 + 6*pow4( msq) + pow4(mU3)) + 16*pow6(m3) + 3*pow6(mU3))
         - pow2(lmQ3MR - logQ3U3)* pow3(m32 - mQ32)*(-((69*mQ32*mU32 +
         pow4(mQ3) - 86*pow4( mU3))*pow5(m3)) + pow3(m3)*(11*mU32*pow4(mQ3) +
         30*mU32*pow4(msq) - 60* msq2*pow4(mU3) + mQ32*(60*msq2*mU32 -
         30*pow4(msq) + 59*pow4(mU3)) - 3* pow6(mQ3) - 99*pow6(mU3)) +
         m3*mU32*(-10*mU32*pow4(mQ3) - 30*mU32*pow4( msq) + mQ32*(-60*msq2*mU32
         + 30*pow4(msq) - 4*pow4(mU3)) + 60*msq2* pow4(mU3) + 3*pow6(mQ3) +
         27*pow6(mU3)) + 18*(mQ32 - mU32)*pow7(m3)) + 2*pow2(mQ32 -
         mU32)*pow2(lmQ3MR - logQ3m3)*pow3(m3)*(-142*m32* mQ32*mU32*(mQ32 +
         mU32) + 87*pow4(mQ3)*pow4(mU3) + pow4(m3)*(220*mQ32* mU32 +
         57*pow4(mQ3) + 57*pow4(mU3)) - 82*(mQ32 + mU32)*pow6(m3) + 27*
         pow8(m3)) - (mQ32 - mU32)*(lmQ3MR - logQ3m3)*(2*(m32 - mQ32)*(m32 -
         mU32)*(mQ32 - mU32)*pow3(m3)*(3*mU32*(10*msq2 + mU32) - 15*m32*(5*
         mQ32 + 4*msq2 + 5*mU32) + mQ32*(30*msq2 + 59*mU32) + 85*pow4(m3) + 3*
         pow4(mQ3)) - 10*(m32 - mQ32)*(m32 - mU32)*(mQ32 - mU32)*(lmQ3MR -
         logQ3mq)*(-11*mQ32*mU32*pow3(m3) + 5*(mQ32 + mU32)*pow5(m3) + pow7(
         m3)) + lmQ3MR*(m32 - mU32)*(68*pow11(m3) + 12*m3*pow4(mU3)*pow6(mQ3) +
         pow5(m3)*(-215*mU32*pow4(mQ3) - 184*mQ32*pow4(mU3) + 3*pow6(mQ3)) +
         pow3(m3)*(105*pow4(mQ3)*pow4(mU3) - 17*mU32*pow6(mQ3)) +
         (353*mQ32*mU32 + 120*pow4(mQ3) + 75*pow4(mU3))*pow7(m3) - (183*mQ32 +
         137*mU32)*pow9( m3)) - (m32 - mQ32)*(lmQ3MR - logQ3U3)*(68*pow11(m3) +
         12*m3* pow4(mQ3)*pow6(mU3) + pow5(m3)*(-184*mU32*pow4(mQ3) -
         215*mQ32*pow4( mU3) + 3*pow6(mU3)) + pow3(m3)*(105*pow4(mQ3)*pow4(mU3)
         - 17*mQ32*pow6( mU3)) + (353*mQ32*mU32 + 75*pow4(mQ3) +
         120*pow4(mU3))*pow7(m3) - (137* mQ32 + 183*mU32)*pow9(m3))) - (m32 -
         mU32)*pow2(lmQ3MR)*(-18*(mQ32 - mU32)*pow11(m3) - 2*((30*msq2 +
         49*mU32)*pow4(mQ3) + 15*mU32*pow4(msq) - mQ32*(30*msq2*mU32 +
         15*pow4(msq) + 58*pow4(mU3)) + 52*pow6(mQ3) - 11*pow6(mU3))*pow7(m3) +
         pow5(m3)*(60*pow4(msq)*pow4(mU3) - 3*pow4(mQ3) *(-20*msq2*mU32 +
         10*pow4(msq) + 9*pow4(mU3)) + (60*msq2 + 169*mU32)* pow6(mQ3) -
         mQ32*(30*mU32*pow4(msq) + 120*msq2*pow4(mU3) + 83*pow6(mU3) ) +
         37*pow8(mQ3)) + (-33*mQ32*mU32 + 86*pow4(mQ3) - 37*pow4(mU3))*pow9(
         m3) + m3*mQ32*(-30*mQ32*msq2*(msq2 + 2*mU32)*pow4(mU3) + 42*pow4(mU3)*
         pow6(mQ3) + pow4(mQ3)*(60*msq2*pow4(mU3) - 19*pow6(mU3)) +
         30*pow4(msq) *pow6(mU3) - 10*mU32*pow8(mQ3) + 3*power10(mQ3)) -
         2*pow3(m3)*(15*mQ32* msq2*(msq2 - 2*mU32)*pow4(mU3) + (60*msq2*mU32 +
         53*pow4(mU3))*pow6( mQ3) + 15*pow4(msq)*pow6(mU3) -
         6*pow4(mQ3)*(5*mU32*pow4(msq) + 5*msq2* pow4(mU3) + 7*pow6(mU3)) +
         17*mU32*pow8(mQ3) + 4*power10(mQ3))) + lmQ3MR*(m32 - mQ32)*(m32 -
         mU32)*(-10*m3*(mQ32 - mU32)*(lmQ3MR - logQ3mq)*pow2(m32 -
         mU32)*(m32*mQ32 + 12*mQ32*msq2 - pow4(mQ3) - 6* pow4(msq)) + 2*(m32 -
         mU32)*(mQ32 - mU32)*(-(mQ32*(61*mQ32 + 30*msq2 + 53*mU32)*pow3(m3)) +
         m3*mQ32*(53*mQ32*mU32 + 30*msq2*mU32 + 3*pow4(mQ3) + 3*pow4(mU3)) +
         (55*mQ32 - 16*mU32)*pow5(m3) + 16*pow7(m3)) - m3*( lmQ3MR -
         logQ3U3)*((-34*mQ32*mU32 + pow4(mQ3) + pow4(mU3))*pow6( m3) +
         38*pow4(mU3)*pow6(mQ3) + pow4(m3)*(29*mU32*pow4(mQ3) + 59*mQ32*
         pow4(mU3) + 9*pow6(mQ3) - pow6(mU3)) + 8*pow4(mQ3)*pow6(mU3) -
         20*mU32* pow8(mQ3) - m32*(88*pow4(mQ3)*pow4(mU3) - 11*mU32*pow6(mQ3) +
         9*mQ32* pow6(mU3) + 10*pow8(mQ3)) +
         6*power10(mQ3))))))/(3.*mQ32*mU32*pow3(m32 - mQ32)*pow3(m32 -
         mU32)*pow5(mQ32 - mU32))+ (1280*Xt6*(-2*mQ32 + 2*mU32 + (mQ32 +
         mU32)*logQ3U3)*pow2(m3*(m32 - mQ32)*mU32*logQ3U3 + (mQ32 -
         mU32)*logQ3m3*pow3(m3)))/( 3.*pow2(m32 - mQ32)*pow2(m32 -
         mU32)*pow5(mQ32 - mU32));
      }

      case Limits::MQ3_EQ_MU3: {
         const double logQ3m3 = std::log(mQ32/m32);
         const double logQ3mq = std::log(mQ32/msq2);
         const double logmqQ3 = -logQ3mq;
         const double logm3Q3 = -logQ3m3;

         return
         (2*Xt4*(120*(m32 - mQ32)*mQ32*(m32 - msq2)*dilog(1 -
         m32/msq2)*(mQ32*msq2 + m32*(-5*mQ32 + 3*msq2) + pow4(m3)) +
         600*msq2*pow4(m3)*pow4(mQ3) + 480*msq2*dilog(1 -
         msq2/mQ32)*pow4(m3)*pow4(mQ3) + 240*msq2*logQ3m3*pow4(m3)*pow4(mQ3) -
         600*msq2*logQ3mq*pow4(m3)*pow4(mQ3) +
         240*msq2*pow2(logQ3mq)*pow4(m3)*pow4(mQ3) + 8*dilog(1 - m32/
         mQ32)*pow2(-(m32*mQ3) + pow3(mQ3))*(-20*m32*mQ32 - 19*pow4(m3) + 15*
         pow4(mQ3)) - 360*mQ32*dilog(1 - msq2/mQ32)*pow4(m3)*pow4(msq) - 180*
         mQ32*pow2(logQ3mq)*pow4(m3)*pow4(msq) + 240*m32*dilog(1 - msq2/
         mQ32)*pow4(mQ3)*pow4(msq) + 120*m32*pow2(logQ3mq)*pow4(mQ3)* pow4(msq)
         - 360*mQ32*msq2*pow6(m3) + 240*mQ32*msq2*dilog(1 - msq2/mQ32)
         *pow6(m3) - 360*mQ32*msq2*logQ3m3*pow6(m3) +
         360*mQ32*msq2*logQ3mq*pow6(m3) + 120*mQ32*msq2*pow2(logQ3mq)*pow6(m3)
         - 10684*pow4(mQ3)*pow6(m3) + 13512*lmQ3MR*pow4(mQ3)*pow6(m3) -
         240*dilog( 1 - msq2/mQ32)*pow4(mQ3)*pow6(m3) +
         19228*logm3Q3*pow4(mQ3)*pow6( m3) +
         11312*lmQ3MR*logm3Q3*pow4(mQ3)*pow6(m3) - 524*logQ3m3
         *pow4(mQ3)*pow6(m3) + 18144*lmQ3MR*logQ3m3*pow4(mQ3)*pow6(m3) +
         24360*logm3Q3*logQ3m3*pow4(mQ3)*pow6(m3) -
         1680*logQ3mq*pow4(mQ3)*pow6(m3) +
         480*lmQ3MR*logQ3mq*pow4(mQ3)*pow6(m3) -
         2020*logm3Q3*logQ3mq*pow4(mQ3)*pow6(m3) -
         1700*logQ3m3*logQ3mq*pow4(mQ3)*pow6(m3) - 4688*pow2(lmQ3MR)*pow4(mQ3)*
         pow6(m3) + 16672*pow2(logQ3m3)*pow4(mQ3)*pow6(m3) -
         120*pow2(logQ3mq)*pow4(mQ3)*pow6(m3) - 120*m32*msq2*pow6(mQ3) -
         720*m32*msq2* dilog(1 - msq2/mQ32)*pow6(mQ3) +
         120*m32*msq2*logQ3m3*pow6(mQ3) + 120*m32*msq2*logQ3mq*pow6(mQ3) -
         360*m32*msq2*pow2(logQ3mq)*pow6(mQ3) + 14474*pow4(m3)*pow6(mQ3) -
         14104*lmQ3MR*pow4(m3)* pow6(mQ3) - 120*dilog(1 -
         msq2/mQ32)*pow4(m3)*pow6(mQ3) - 13712*logm3Q3*pow4(m3)*pow6(mQ3) -
         8064*lmQ3MR*logm3Q3*pow4(m3)*pow6( mQ3) +
         2504*logQ3m3*pow4(m3)*pow6(mQ3) -
         11240*lmQ3MR*logQ3m3*pow4(m3)*pow6(mQ3) +
         1728*logm3Q3*logQ3m3*pow4(m3)* pow6(mQ3) +
         2000*logQ3mq*pow4(m3)*pow6(mQ3) -
         320*lmQ3MR*logQ3mq*pow4(m3)*pow6(mQ3) + 1440*logm3Q3*logQ3mq*pow4(
         m3)*pow6(mQ3) + 1000*logQ3m3*logQ3mq*pow4(m3)*pow6(mQ3) +
         4792*pow2(lmQ3MR)*pow4(m3)*pow6(mQ3) + 2548*pow2(logQ3m3)*pow4(
         m3)*pow6(mQ3) - 60*pow2(logQ3mq)*pow4(m3)*pow6(mQ3) + 120*dilog( 1 -
         msq2/mQ32)*pow4(msq)*pow6(mQ3) + 60*pow2(logQ3mq)*pow4(msq)* pow6(mQ3)
         + 2463*mQ32*pow8(m3) - 5440*lmQ3MR*mQ32*pow8(m3) - 5680*mQ32*
         logm3Q3*pow8(m3) - 7168*lmQ3MR*mQ32*logm3Q3*pow8(m3) - 248*
         mQ32*logQ3m3*pow8(m3) - 11720*lmQ3MR*mQ32*logQ3m3*pow8(m3) -
         21312*mQ32*logm3Q3*logQ3m3*pow8(m3) + 640*mQ32*logQ3mq*pow8(m3) -
         320*lmQ3MR*mQ32*logQ3mq*pow8(m3) + 1280*mQ32* logm3Q3*logQ3mq*pow8(m3)
         + 1480*mQ32*logQ3m3*logQ3mq*pow8(m3) + 2292*mQ32*pow2(lmQ3MR)*pow8(m3)
         - 18996*mQ32* pow2(logQ3m3)*pow8(m3) - 5868*m32*pow8(mQ3) +
         6312*lmQ3MR*m32* pow8(mQ3) - 120*msq2*pow8(mQ3) + 480*m32*dilog(1 -
         msq2/mQ32)*pow8(mQ3) - 630*m32*logm3Q3*pow8(mQ3) +
         2184*lmQ3MR*m32*logm3Q3*pow8( mQ3) - 2142*m32*logQ3m3*pow8(mQ3) +
         2280*lmQ3MR*m32*logQ3m3 *pow8(mQ3) - 468*m32*logm3Q3*logQ3m3*pow8(mQ3)
         - 1120*m32* logQ3mq*pow8(mQ3) + 80*lmQ3MR*m32*logQ3mq*pow8(mQ3) +
         120*msq2*logQ3mq*pow8(mQ3) - 390*m32*logm3Q3*logQ3mq*pow8(mQ3) -
         390*m32*logQ3m3*logQ3mq*pow8(mQ3) - 2448*m32*pow2(lmQ3MR)*pow8(mQ3) -
         180*m32*pow2(logQ3m3)*pow8(mQ3) + 240*m32*pow2(logQ3mq)*pow8(mQ3) -
         176*power10(m3) + 624* lmQ3MR*power10(m3) + 794*logm3Q3*power10(m3) +
         1736*lmQ3MR*logm3Q3*power10(m3) + 266*logQ3m3*power10(m3) +
         2536*lmQ3MR* logQ3m3*power10(m3) + 6572*logm3Q3*logQ3m3*power10( m3) -
         80*logQ3mq*power10(m3) + 80*lmQ3MR*logQ3mq*power10( m3) -
         310*logm3Q3*logQ3mq*power10(m3) - 390*logQ3m3* logQ3mq*power10(m3) -
         448*pow2(lmQ3MR)*power10(m3) + 6220*pow2( logQ3m3)*power10(m3) -
         209*power10(mQ3) - 904*lmQ3MR*power10(mQ3) - 120*dilog(1 -
         msq2/mQ32)*power10(mQ3) + 144*logQ3m3*power10( mQ3) +
         240*logQ3mq*power10(mQ3) + 500*pow2(lmQ3MR)*power10(mQ3) -
         72*pow2(logQ3m3)*power10(mQ3) - 60*pow2(logQ3mq)*
         power10(mQ3)))/(9.*pow4(m32 - mQ32)*pow6(mQ3))+
         (-32*m3*Xt5*(60*mQ32*(-m32 + mQ32)*(m32 - msq2)*(mQ32 - msq2)*dilog(1
         - m32/msq2) + 2*mQ32*dilog(1 - m32/mQ32)*pow3(m32 - mQ32) -
         60*mQ32*msq2* pow4(m3) + 60*mQ32*msq2*dilog(1 - msq2/mQ32)*pow4(m3) -
         60*mQ32*msq2* logQ3m3*pow4(m3) + 90*mQ32*msq2*logQ3mq*pow4(m3) + 30*
         mQ32*msq2*logmqQ3*pow4(m3) + 30*mQ32*msq2*pow2(logQ3mq)* pow4(m3) +
         120*m32*msq2*pow4(mQ3) + 60*m32*msq2*logQ3m3*pow4(mQ3) -
         180*m32*msq2*logQ3mq*pow4(mQ3) - 60*m32*msq2*logmqQ3* pow4(mQ3) -
         345*pow4(m3)*pow4(mQ3) + 274*lmQ3MR*pow4(m3)*pow4(mQ3) - 60*dilog(1 -
         msq2/mQ32)*pow4(m3)*pow4(mQ3) - 342*logm3Q3*pow4(m3) *pow4(mQ3) -
         538*logQ3m3*pow4(m3)*pow4(mQ3) + 220*lmQ3MR*logQ3m3*pow4(m3)*pow4(mQ3)
         - 492*logm3Q3*logQ3m3*pow4(m3) *pow4(mQ3) -
         15*logQ3mq*pow4(m3)*pow4(mQ3) + 50*logQ3m3* logQ3mq*pow4(m3)*pow4(mQ3)
         - 602*pow2(logQ3m3)*pow4(m3)* pow4(mQ3) -
         30*pow2(logQ3mq)*pow4(m3)*pow4(mQ3) - 60*m32*mQ32* dilog(1 -
         msq2/mQ32)*pow4(msq) + 60*m32*mQ32*logQ3mq*pow4(msq) +
         60*m32*mQ32*logmqQ3*pow4(msq) - 30*m32*mQ32*pow2(logQ3mq) *pow4(msq) -
         30*logQ3mq*pow4(m3)*pow4(msq) - 30*logmqQ3* pow4(m3)*pow4(msq) +
         60*dilog(1 - msq2/mQ32)*pow4(mQ3)*pow4(msq) - 30*
         logQ3mq*pow4(mQ3)*pow4(msq) - 30*logmqQ3*pow4(mQ3)*pow4( msq) +
         30*pow2(logQ3mq)*pow4(mQ3)*pow4(msq) + 163*mQ32*pow6(m3) -
         122*lmQ3MR*mQ32*pow6(m3) + 70*mQ32*logm3Q3*pow6(m3) + 277*mQ32*
         logQ3m3*pow6(m3) - 156*lmQ3MR*mQ32*logQ3m3*pow6(m3) + 408*
         mQ32*logm3Q3*logQ3m3*pow6(m3) + 5*mQ32*logQ3mq*pow6( m3) +
         5*mQ32*logQ3m3*logQ3mq*pow6(m3) + 495*mQ32*pow2(logQ3m3)*pow6(m3) +
         265*m32*pow6(mQ3) - 246*lmQ3MR*m32*pow6(mQ3) - 60*msq2*pow6(mQ3) +
         60*m32*dilog(1 - msq2/mQ32)*pow6(mQ3) - 60*msq2* dilog(1 -
         msq2/mQ32)*pow6(mQ3) + 274*m32*logm3Q3*pow6(mQ3) + 317*
         m32*logQ3m3*pow6(mQ3) - 96*lmQ3MR*m32*logQ3m3*pow6(mQ3) +
         15*m32*logQ3mq*pow6(mQ3) + 90*msq2*logQ3mq*pow6(mQ3) -
         55*m32*logQ3m3*logQ3mq*pow6(mQ3) + 30*msq2*logmqQ3* pow6(mQ3) +
         87*m32*pow2(logQ3m3)*pow6(mQ3) + 30*m32*pow2(logQ3mq)*pow6(mQ3) -
         30*msq2*pow2(logQ3mq)*pow6(mQ3) - 16* pow8(m3) + 16*lmQ3MR*pow8(m3) -
         2*logm3Q3*pow8(m3) - 50*logQ3m3*pow8(m3) + 32*lmQ3MR*logQ3m3*pow8(m3)
         - 124*logm3Q3* logQ3m3*pow8(m3) - 156*pow2(logQ3m3)*pow8(m3) -
         67*pow8( mQ3) + 78*lmQ3MR*pow8(mQ3) - 6*logQ3m3*pow8(mQ3) -
         5*logQ3mq*pow8(mQ3)))/(9.*pow4(m32 - mQ32)*pow6(mQ3))+
         (640*Xt6*pow2(-(m3*mQ32) + pow3(m3) + logQ3m3*pow3(m3)))/(9.*pow4(-
         (m32*mQ3) + pow3(mQ3)));
      }

      case Limits::MQ3_EQ_M3: {
         const double logQ3U3 = std::log(mQ32/mU32);
         const double logQ3mq = std::log(mQ32/msq2);

         return
         -(Xt4*(122320*msq2*mU32*pow14(mQ3) - 6560*lmQ3MR*msq2*mU32*pow14(mQ3)
         + 40992*msq2*mU32*dilog(1 - mQ32/mU32)*pow14(mQ3) - 5376*lmQ3MR*msq2*
         mU32*dilog(1 - mQ32/mU32)*pow14(mQ3) + 5760*msq2*mU32*dilog(1 - msq2/
         mU32)*pow14(mQ3) - 5280*msq2*mU32*logQ3mq*pow14(mQ3) - 3840*
         lmQ3MR*msq2*mU32*logQ3mq*pow14(mQ3) + 960*msq2*mU32*dilog(1 -
         mQ32/mU32)*logQ3mq*pow14(mQ3) - 45464*msq2*mU32*logQ3U3* pow14(mQ3) -
         26288*lmQ3MR*msq2*mU32*logQ3U3*pow14(mQ3) - 8208* msq2*mU32*dilog(1 -
         mQ32/mU32)*logQ3U3*pow14(mQ3) - 2880*msq2* mU32*dilog(1 -
         msq2/mU32)*logQ3U3*pow14(mQ3) + 720*msq2*mU32*
         logQ3mq*logQ3U3*pow14(mQ3) -
         1920*lmQ3MR*msq2*mU32*logQ3mq*logQ3U3*pow14(mQ3) + 192*msq2*pow16(mQ3)
         + 5184* lmQ3MR*msq2*pow16(mQ3) - 130960*mU32*pow16(mQ3) +
         6560*lmQ3MR*mU32* pow16(mQ3) - 40992*mU32*dilog(1 -
         mQ32/mU32)*pow16(mQ3) + 5376*lmQ3MR* mU32*dilog(1 -
         mQ32/mU32)*pow16(mQ3) - 960*msq2*logQ3mq*pow16( mQ3) +
         960*lmQ3MR*msq2*logQ3mq*pow16(mQ3) + 12000*mU32*logQ3mq*pow16(mQ3) +
         3840*lmQ3MR*mU32*logQ3mq*pow16(mQ3) - 960* mU32*dilog(1 -
         mQ32/mU32)*logQ3mq*pow16(mQ3) - 96*msq2*logQ3U3*pow16(mQ3) +
         96*lmQ3MR*msq2*logQ3U3*pow16(mQ3) + 49784* mU32*logQ3U3*pow16(mQ3) +
         26288*lmQ3MR*mU32*logQ3U3* pow16(mQ3) + 8208*mU32*dilog(1 -
         mQ32/mU32)*logQ3U3*pow16(mQ3) - 9840*mU32*logQ3mq*logQ3U3*pow16(mQ3) +
         1920*lmQ3MR*mU32* logQ3mq*logQ3U3*pow16(mQ3) - 192*pow18(mQ3) - 5184*
         lmQ3MR*pow18(mQ3) + 960*logQ3mq*pow18(mQ3) -
         960*lmQ3MR*logQ3mq*pow18(mQ3) + 96*logQ3U3*pow18(mQ3) -
         96*lmQ3MR*logQ3U3*pow18(mQ3) - 2496*msq2*mU32*pow14(mQ3)*pow2(lmQ3MR)
         + 22752* msq2*mU32*logQ3U3*pow14(mQ3)*pow2(lmQ3MR) -
         5376*msq2*pow16(mQ3) *pow2(lmQ3MR) + 2496*mU32*pow16(mQ3)*pow2(lmQ3MR)
         - 22752*mU32*logQ3U3*pow16(mQ3)*pow2(lmQ3MR) +
         5376*pow18(mQ3)*pow2(lmQ3MR) + 3330*msq2*mU32*pow14(mQ3)*pow2(logQ3mq)
         - 1665*msq2*mU32*logQ3U3*pow14(mQ3)*pow2(logQ3mq) -
         270*mU32*pow16(mQ3)*pow2( logQ3mq) +
         135*mU32*logQ3U3*pow16(mQ3)*pow2(logQ3mq) +
         21480*msq2*mU32*pow14(mQ3)*pow2(logQ3U3) - 2496*lmQ3MR*
         msq2*mU32*pow14(mQ3)*pow2(logQ3U3) +
         3840*msq2*mU32*logQ3mq*pow14(mQ3)*pow2(logQ3U3) -
         18600*mU32*pow16(mQ3)*pow2(logQ3U3) +
         2496*lmQ3MR*mU32*pow16(mQ3)*pow2(logQ3U3) - 960*
         mU32*logQ3mq*pow16(mQ3)*pow2(logQ3U3) - 5640*msq2*mU32*
         pow14(mQ3)*pow3(logQ3U3) + 4200*mU32*pow16(mQ3)*pow3(logQ3U3) -
         29968*msq2*pow12(mU3)*pow4(mQ3) - 58720*lmQ3MR*msq2*pow12(mU3)
         *pow4(mQ3) + 3936*msq2*dilog(1 - mQ32/mU32)*pow12(mU3)*pow4(mQ3) +
         5376*lmQ3MR*msq2*dilog(1 - mQ32/mU32)*pow12(mU3)*pow4(mQ3) -
         2880*msq2* dilog(1 - msq2/mU32)*pow12(mU3)*pow4(mQ3) +
         17760*msq2*logQ3mq* pow12(mU3)*pow4(mQ3) -
         960*lmQ3MR*msq2*logQ3mq*pow12(mU3)*pow4( mQ3) - 960*msq2*dilog(1 -
         mQ32/mU32)*logQ3mq*pow12(mU3)*pow4( mQ3) +
         520*msq2*logQ3U3*pow12(mU3)*pow4(mQ3) - 63920*lmQ3MR*
         msq2*logQ3U3*pow12(mU3)*pow4(mQ3) - 528*msq2*dilog(1 - mQ32/
         mU32)*logQ3U3*pow12(mU3)*pow4(mQ3) - 1440*msq2*dilog(1 - msq2/
         mU32)*logQ3U3*pow12(mU3)*pow4(mQ3) + 13200*msq2*logQ3mq*
         logQ3U3*pow12(mU3)*pow4(mQ3) + 29376*msq2*pow12(mU3)*pow2(
         lmQ3MR)*pow4(mQ3) + 12000*msq2*logQ3U3*pow12(mU3)*pow2(lmQ3MR)*
         pow4(mQ3) - 1890*msq2*pow12(mU3)*pow2(logQ3mq)*pow4(mQ3) - 945*
         msq2*logQ3U3*pow12(mU3)*pow2(logQ3mq)*pow4(mQ3) + 14616*
         msq2*pow12(mU3)*pow2(logQ3U3)*pow4(mQ3) - 17376*lmQ3MR*msq2*
         pow12(mU3)*pow2(logQ3U3)*pow4(mQ3) + 2400*msq2*logQ3mq*
         pow12(mU3)*pow2(logQ3U3)*pow4(mQ3) + 4056*msq2*pow12(mU3)*pow3(
         logQ3U3)*pow4(mQ3) + 8640*mU32*pow12(mQ3)*pow4(msq) - 14400*
         mU32*dilog(1 - msq2/mU32)*pow12(mQ3)*pow4(msq) -
         8640*mU32*logQ3mq*pow12(mQ3)*pow4(msq) -
         4320*mU32*logQ3U3*pow12(mQ3)*pow4( msq) + 7200*mU32*dilog(1 -
         msq2/mU32)*logQ3U3*pow12(mQ3)*pow4( msq) +
         18720*mU32*logQ3mq*logQ3U3*pow12(mQ3)*pow4(msq) -
         7290*mU32*pow12(mQ3)*pow2(logQ3mq)*pow4(msq) +
         3645*mU32*logQ3U3*pow12(mQ3)*pow2(logQ3mq)*pow4(msq) + 90*mQ32*pow12(
         mU3)*pow2(logQ3mq)*pow4(msq) + 45*mQ32*logQ3U3*pow12(mU3)
         *pow2(logQ3mq)*pow4(msq) -
         7200*mU32*pow12(mQ3)*pow2(logQ3U3)*pow4(msq) -
         7200*mU32*logQ3mq*pow12(mQ3)*pow2(logQ3U3)*pow4(msq) +
         3600*mU32*pow12(mQ3)*pow3(logQ3U3)*pow4(msq) -
         536144*msq2*pow12(mQ3)*pow4(mU3) - 85088*lmQ3MR*msq2*pow12(mQ3)*pow4(
         mU3) - 155424*msq2*dilog(1 - mQ32/mU32)*pow12(mQ3)*pow4(mU3) + 26880*
         lmQ3MR*msq2*dilog(1 - mQ32/mU32)*pow12(mQ3)*pow4(mU3) + 11520*msq2*
         dilog(1 - msq2/mU32)*pow12(mQ3)*pow4(mU3) + 60000*msq2*logQ3mq*
         pow12(mQ3)*pow4(mU3) + 4800*lmQ3MR*msq2*logQ3mq*pow12(mQ3)*pow4( mU3)
         - 4800*msq2*dilog(1 - mQ32/mU32)*logQ3mq*pow12(mQ3)*pow4( mU3) +
         13032*msq2*logQ3U3*pow12(mQ3)*pow4(mU3) + 119760*lmQ3MR*
         msq2*logQ3U3*pow12(mQ3)*pow4(mU3) - 14928*msq2*dilog(1 - mQ32/
         mU32)*logQ3U3*pow12(mQ3)*pow4(mU3) - 11520*msq2*dilog(1 - msq2/
         mU32)*logQ3U3*pow12(mQ3)*pow4(mU3) - 30960*msq2*logQ3mq*
         logQ3U3*pow12(mQ3)*pow4(mU3) + 7680*lmQ3MR*msq2*logQ3mq*
         logQ3U3*pow12(mQ3)*pow4(mU3) + 559184*pow14(mQ3)*pow4(mU3) +
         85088*lmQ3MR*pow14(mQ3)*pow4(mU3) + 155424*dilog(1 - mQ32/mU32)*pow14(
         mQ3)*pow4(mU3) - 26880*lmQ3MR*dilog(1 -
         mQ32/mU32)*pow14(mQ3)*pow4(mU3) - 5760*dilog(1 -
         msq2/mU32)*pow14(mQ3)*pow4(mU3) - 73440*logQ3mq *pow14(mQ3)*pow4(mU3)
         - 4800*lmQ3MR*logQ3mq*pow14(mQ3)*pow4(mU3) + 4800*dilog(1 -
         mQ32/mU32)*logQ3mq*pow14(mQ3)*pow4(mU3) - 7272*
         logQ3U3*pow14(mQ3)*pow4(mU3) - 119760*lmQ3MR*logQ3U3*
         pow14(mQ3)*pow4(mU3) + 14928*dilog(1 - mQ32/mU32)*logQ3U3*pow14(
         mQ3)*pow4(mU3) + 2880*dilog(1 - msq2/mU32)*logQ3U3*pow14(mQ3)*
         pow4(mU3) + 25200*logQ3mq*logQ3U3*pow14(mQ3)*pow4(mU3) -
         7680*lmQ3MR*logQ3mq*logQ3U3*pow14(mQ3)*pow4(mU3) + 93120*
         msq2*pow12(mQ3)*pow2(lmQ3MR)*pow4(mU3) - 79008*msq2*logQ3U3*
         pow12(mQ3)*pow2(lmQ3MR)*pow4(mU3) -
         93120*pow14(mQ3)*pow2(lmQ3MR)*pow4( mU3) +
         79008*logQ3U3*pow14(mQ3)*pow2(lmQ3MR)*pow4(mU3) + 3510*
         msq2*pow12(mQ3)*pow2(logQ3mq)*pow4(mU3) -
         5085*msq2*logQ3U3*pow12(mQ3)*pow2(logQ3mq)*pow4(mU3) -
         1530*pow14(mQ3)*pow2( logQ3mq)*pow4(mU3) +
         1035*logQ3U3*pow14(mQ3)*pow2(logQ3mq)*pow4(mU3) -
         27816*msq2*pow12(mQ3)*pow2(logQ3U3)* pow4(mU3) +
         18240*lmQ3MR*msq2*pow12(mQ3)*pow2(logQ3U3)*pow4(mU3) +
         1920*msq2*logQ3mq*pow12(mQ3)*pow2(logQ3U3)*pow4(mU3) +
         26376*pow14(mQ3)*pow2(logQ3U3)*pow4(mU3) - 18240*lmQ3MR*pow14(
         mQ3)*pow2(logQ3U3)*pow4(mU3) + 6720*logQ3mq*pow14(mQ3)*
         pow2(logQ3U3)*pow4(mU3) -
         16200*msq2*pow12(mQ3)*pow3(logQ3U3)*pow4(mU3) +
         11880*pow14(mQ3)*pow3(logQ3U3)*pow4(mU3) + 90*(3*mQ32 +
         msq2)*mU32*dilog(1 - msq2/mQ32)*(-2*mQ32 + 2*mU32 + (mQ32 +
         mU32)*logQ3U3)*pow2(mQ32 - msq2)*pow4(mQ32 - mU32) + 29968*
         pow12(mU3)*pow6(mQ3) + 58720*lmQ3MR*pow12(mU3)*pow6(mQ3) -
         3936*dilog(1 - mQ32/mU32)*pow12(mU3)*pow6(mQ3) - 5376*lmQ3MR*dilog(1 -
         mQ32/mU32)* pow12(mU3)*pow6(mQ3) + 2880*dilog(1 -
         msq2/mU32)*pow12(mU3)*pow6(mQ3) - 15840*logQ3mq*pow12(mU3)*pow6(mQ3) +
         960*lmQ3MR*logQ3mq* pow12(mU3)*pow6(mQ3) + 960*dilog(1 -
         mQ32/mU32)*logQ3mq*pow12( mU3)*pow6(mQ3) -
         520*logQ3U3*pow12(mU3)*pow6(mQ3) + 63920*
         lmQ3MR*logQ3U3*pow12(mU3)*pow6(mQ3) + 528*dilog(1 - mQ32/mU32)*
         logQ3U3*pow12(mU3)*pow6(mQ3) + 1440*dilog(1 -
         msq2/mU32)*logQ3U3*pow12(mU3)*pow6(mQ3) - 12240*logQ3mq*logQ3U3*
         pow12(mU3)*pow6(mQ3) - 29376*pow12(mU3)*pow2(lmQ3MR)*pow6(mQ3) -
         12000* logQ3U3*pow12(mU3)*pow2(lmQ3MR)*pow6(mQ3) + 1710*pow12(mU3)*
         pow2(logQ3mq)*pow6(mQ3) + 855*logQ3U3*pow12(mU3)*pow2(
         logQ3mq)*pow6(mQ3) - 14616*pow12(mU3)*pow2(logQ3U3)*pow6( mQ3) +
         17376*lmQ3MR*pow12(mU3)*pow2(logQ3U3)*pow6(mQ3) - 2400*
         logQ3mq*pow12(mU3)*pow2(logQ3U3)*pow6(mQ3) - 4056*pow12(
         mU3)*pow3(logQ3U3)*pow6(mQ3) + 90*pow12(mU3)*pow2(logQ3mq )*pow6(msq)
         + 45*logQ3U3*pow12(mU3)*pow2(logQ3mq)*pow6( msq) -
         919840*pow12(mQ3)*pow6(mU3) - 301504*lmQ3MR*pow12(mQ3)*pow6(mU3) -
         216384*dilog(1 - mQ32/mU32)*pow12(mQ3)*pow6(mU3) + 53760*lmQ3MR*
         dilog(1 - mQ32/mU32)*pow12(mQ3)*pow6(mU3) + 2880*dilog(1 - msq2/mU32)*
         pow12(mQ3)*pow6(mU3) + 154560*logQ3mq*pow12(mQ3)*pow6(mU3) -
         9600*dilog(1 - mQ32/mU32)*logQ3mq*pow12(mQ3)*pow6(mU3) - 284624*
         logQ3U3*pow12(mQ3)*pow6(mU3) + 137056*lmQ3MR*logQ3U3*
         pow12(mQ3)*pow6(mU3) - 89952*dilog(1 - mQ32/mU32)*logQ3U3*pow12(
         mQ3)*pow6(mU3) + 4320*dilog(1 - msq2/mU32)*logQ3U3*pow12(mQ3)*
         pow6(mU3) - 4320*logQ3mq*logQ3U3*pow12(mQ3)*pow6(mU3) +
         11520*lmQ3MR*logQ3mq*logQ3U3*pow12(mQ3)*pow6(mU3) +
         240000*pow12(mQ3)*pow2(lmQ3MR)*pow6(mU3) - 88512*logQ3U3*pow12(
         mQ3)*pow2(lmQ3MR)*pow6(mU3) - 1260*pow12(mQ3)*pow2(logQ3mq)* pow6(mU3)
         + 2430*logQ3U3*pow12(mQ3)*pow2(logQ3mq)*pow6( mU3) -
         41424*pow12(mQ3)*pow2(logQ3U3)*pow6(mU3) + 62496*lmQ3MR*
         pow12(mQ3)*pow2(logQ3U3)*pow6(mU3) - 11040*logQ3mq*pow12(
         mQ3)*pow2(logQ3U3)*pow6(mU3) -
         39888*pow12(mQ3)*pow3(logQ3U3)*pow6(mU3) + 2880*dilog(1 -
         msq2/mU32)*pow6(mQ3)*pow6(msq)*pow6( mU3) + 4320*dilog(1 -
         msq2/mU32)*logQ3U3*pow6(mQ3)*pow6(msq)* pow6(mU3) -
         2880*logQ3mq*logQ3U3*pow6(mQ3)*pow6(msq)* pow6(mU3) +
         540*pow2(logQ3mq)*pow6(mQ3)*pow6(msq)*pow6(mU3) +
         2250*logQ3U3*pow2(logQ3mq)*pow6(mQ3)*pow6(msq)*pow6(mU3) +
         1440*pow2(logQ3U3)*pow6(mQ3)*pow6(msq)*pow6(mU3) -
         4320*logQ3mq*pow2(logQ3U3)*pow6(mQ3)*pow6(msq)*pow6(mU3) + 2160*
         pow3(logQ3U3)*pow6(mQ3)*pow6(msq)*pow6(mU3) - 14400*dilog(1 -
         msq2/mU32)*pow4(mU3)*pow6(msq)*pow8(mQ3) - 1440*dilog(1 - msq2/mU32)*
         logQ3U3*pow4(mU3)*pow6(msq)*pow8(mQ3) + 14400*logQ3mq*
         logQ3U3*pow4(mU3)*pow6(msq)*pow8(mQ3) - 6750*pow2(logQ3mq
         )*pow4(mU3)*pow6(msq)*pow8(mQ3) -
         855*logQ3U3*pow2(logQ3mq)*pow4(mU3)*pow6(msq)*pow8(mQ3) -
         7200*pow2(logQ3U3)*pow4( mU3)*pow6(msq)*pow8(mQ3) +
         1440*logQ3mq*pow2(logQ3U3)* pow4(mU3)*pow6(msq)*pow8(mQ3) -
         720*pow3(logQ3U3)*pow4(mU3)* pow6(msq)*pow8(mQ3) +
         17280*pow4(msq)*pow6(mU3)*pow8(mQ3) + 25920* dilog(1 -
         msq2/mU32)*pow4(msq)*pow6(mU3)*pow8(mQ3) -
         17280*logQ3mq*pow4(msq)*pow6(mU3)*pow8(mQ3) + 20160*logQ3U3*pow4(msq)*
         pow6(mU3)*pow8(mQ3) - 7200*dilog(1 - msq2/mU32)*logQ3U3*pow4(
         msq)*pow6(mU3)*pow8(mQ3) - 31680*logQ3mq*logQ3U3*pow4(
         msq)*pow6(mU3)*pow8(mQ3) + 12060*pow2(logQ3mq)*pow4(msq)*pow6(
         mU3)*pow8(mQ3) - 3510*logQ3U3*pow2(logQ3mq)*pow4(msq)*
         pow6(mU3)*pow8(mQ3) + 14400*pow2(logQ3U3)*pow4(msq)*pow6(mU3)*
         pow8(mQ3) + 7200*logQ3mq*pow2(logQ3U3)*pow4(msq)*pow6( mU3)*pow8(mQ3)
         - 3600*pow3(logQ3U3)*pow4(msq)*pow6(mU3)*pow8( mQ3) - 20160*dilog(1 -
         msq2/mU32)*pow4(msq)*pow6(mQ3)*pow8(mU3) - 5760*
         logQ3U3*pow4(msq)*pow6(mQ3)*pow8(mU3) - 10080*dilog(1 - msq2/
         mU32)*logQ3U3*pow4(msq)*pow6(mQ3)*pow8(mU3) +
         23040*logQ3mq*logQ3U3*pow4(msq)*pow6(mQ3)*pow8(mU3) -
         9180*pow2(logQ3mq)*pow4(msq)*pow6(mQ3)*pow8(mU3) - 4950*logQ3U3*pow2(
         logQ3mq)*pow4(msq)*pow6(mQ3)*pow8(mU3) -
         14400*pow2(logQ3U3)*pow4(msq)*pow6(mQ3)*pow8(mU3) +
         10080*logQ3mq*pow2(logQ3U3)*pow4(msq)*pow6(mQ3)*pow8(mU3) -
         5040*pow3(logQ3U3)* pow4(msq)*pow6(mQ3)*pow8(mU3) + 2880*dilog(1 -
         msq2/mU32)*pow4(mQ3)* pow6(msq)*pow8(mU3) + 1440*dilog(1 -
         msq2/mU32)*logQ3U3*pow4( mQ3)*pow6(msq)*pow8(mU3) -
         2880*logQ3mq*logQ3U3*pow4(mQ3) *pow6(msq)*pow8(mU3) +
         2340*pow2(logQ3mq)*pow4(mQ3)*pow6(msq)* pow8(mU3) +
         810*logQ3U3*pow2(logQ3mq)*pow4(mQ3)*pow6(msq) *pow8(mU3) +
         1440*pow2(logQ3U3)*pow4(mQ3)*pow6(msq)*pow8(mU3) -
         1440*logQ3mq*pow2(logQ3U3)*pow4(mQ3)*pow6(msq)*pow8(mU3) +
         720*pow3(logQ3U3)*pow4(mQ3)*pow6(msq)*pow8(mU3) - 716896*msq2*
         pow8(mQ3)*pow8(mU3) - 406144*lmQ3MR*msq2*pow8(mQ3)*pow8(mU3) - 126528*
         msq2*dilog(1 - mQ32/mU32)*pow8(mQ3)*pow8(mU3) + 53760*lmQ3MR*msq2*
         dilog(1 - mQ32/mU32)*pow8(mQ3)*pow8(mU3) + 2880*msq2*dilog(1 - msq2/
         mU32)*pow8(mQ3)*pow8(mU3) + 176640*msq2*logQ3mq*pow8(mQ3)*pow8( mU3) -
         4800*lmQ3MR*msq2*logQ3mq*pow8(mQ3)*pow8(mU3) - 9600*msq2* dilog(1 -
         mQ32/mU32)*logQ3mq*pow8(mQ3)*pow8(mU3) - 385168*msq2*
         logQ3U3*pow8(mQ3)*pow8(mU3) - 44416*lmQ3MR*msq2*logQ3U3*
         pow8(mQ3)*pow8(mU3) - 102816*msq2*dilog(1 - mQ32/mU32)*logQ3U3*
         pow8(mQ3)*pow8(mU3) + 12960*msq2*dilog(1 - msq2/mU32)*logQ3U3*
         pow8(mQ3)*pow8(mU3) + 21600*msq2*logQ3mq*logQ3U3*pow8( mQ3)*pow8(mU3)
         + 7680*lmQ3MR*msq2*logQ3mq*logQ3U3*pow8( mQ3)*pow8(mU3) +
         266880*msq2*pow2(lmQ3MR)*pow8(mQ3)*pow8(mU3) - 19008*
         msq2*logQ3U3*pow2(lmQ3MR)*pow8(mQ3)*pow8(mU3) - 3060*msq2*pow2(
         logQ3mq)*pow8(mQ3)*pow8(mU3) + 6030*msq2*logQ3U3*pow2(
         logQ3mq)*pow8(mQ3)*pow8(mU3) - 54768*msq2*pow2(logQ3U3)*
         pow8(mQ3)*pow8(mU3) + 62112*lmQ3MR*msq2*pow2(logQ3U3)*pow8(mQ3)*
         pow8(mU3) - 12000*msq2*logQ3mq*pow2(logQ3U3)*pow8(mQ3)* pow8(mU3) -
         16080*msq2*pow3(logQ3U3)*pow8(mQ3)*pow8(mU3) -
         23040*pow4(msq)*pow4(mU3)*power10(mQ3) + 8640*dilog(1 - msq2/mU32)*
         pow4(msq)*pow4(mU3)*power10(mQ3) + 23040*logQ3mq*pow4(msq)*pow4(
         mU3)*power10(mQ3) - 5760*logQ3U3*pow4(msq)*pow4(mU3)*power10( mQ3) +
         10080*dilog(1 - msq2/mU32)*logQ3U3*pow4(msq)*pow4(mU3)* power10(mQ3) -
         11520*logQ3mq*logQ3U3*pow4(msq)*pow4(mU3)* power10(mQ3) +
         4770*pow2(logQ3mq)*pow4(msq)*pow4(mU3)*power10( mQ3) +
         4905*logQ3U3*pow2(logQ3mq)*pow4(msq)*pow4(mU3)* power10(mQ3) +
         8640*pow2(logQ3U3)*pow4(msq)*pow4(mU3)*power10( mQ3) -
         10080*logQ3mq*pow2(logQ3U3)*pow4(msq)*pow4(mU3)* power10(mQ3) +
         5040*pow3(logQ3U3)*pow4(msq)*pow4(mU3)*power10( mQ3) +
         8640*mU32*dilog(1 - msq2/mU32)*pow6(msq)*power10(mQ3) - 4320*
         mU32*dilog(1 - msq2/mU32)*logQ3U3*pow6(msq)*power10(mQ3) - 8640*
         mU32*logQ3mq*logQ3U3*pow6(msq)*power10(mQ3) + 4230*mU32*
         pow2(logQ3mq)*pow6(msq)*power10(mQ3) - 2115*mU32*logQ3U3*
         pow2(logQ3mq)*pow6(msq)*power10(mQ3) +
         4320*mU32*pow2(logQ3U3)*pow6(msq)*power10(mQ3) +
         4320*mU32*logQ3mq*pow2(logQ3U3)*pow6(msq)*power10(mQ3) -
         2160*mU32*pow3(logQ3U3)*pow6( msq)*power10(mQ3) +
         902560*msq2*pow6(mU3)*power10(mQ3) + 301504*lmQ3MR*
         msq2*pow6(mU3)*power10(mQ3) + 216384*msq2*dilog(1 - mQ32/mU32)*pow6(
         mU3)*power10(mQ3) - 53760*lmQ3MR*msq2*dilog(1 - mQ32/mU32)*pow6(mU3)*
         power10(mQ3) - 31680*msq2*dilog(1 - msq2/mU32)*pow6(mU3)*power10(mQ3)
         - 156480*msq2*logQ3mq*pow6(mU3)*power10(mQ3) + 9600*msq2*dilog(1 -
         mQ32/mU32)*logQ3mq*pow6(mU3)*power10(mQ3) +
         264464*msq2*logQ3U3*pow6(mU3)*power10(mQ3) -
         137056*lmQ3MR*msq2*logQ3U3* pow6(mU3)*power10(mQ3) +
         89952*msq2*dilog(1 - mQ32/mU32)*logQ3U3 *pow6(mU3)*power10(mQ3) -
         1440*msq2*dilog(1 - msq2/mU32)*logQ3U3 *pow6(mU3)*power10(mQ3) +
         40800*msq2*logQ3mq*logQ3U3* pow6(mU3)*power10(mQ3) -
         11520*lmQ3MR*msq2*logQ3mq*logQ3U3*pow6(mU3)*power10(mQ3) -
         240000*msq2*pow2(lmQ3MR)*pow6(mU3)* power10(mQ3) +
         88512*msq2*logQ3U3*pow2(lmQ3MR)*pow6(mU3)* power10(mQ3) -
         11340*msq2*pow2(logQ3mq)*pow6(mU3)*power10(mQ3) -
         1170*msq2*logQ3U3*pow2(logQ3mq)*pow6(mU3)*power10(mQ3) +
         25584*msq2*pow2(logQ3U3)*pow6(mU3)*power10(mQ3) - 62496*lmQ3MR*
         msq2*pow2(logQ3U3)*pow6(mU3)*power10(mQ3) +
         8160*msq2*logQ3mq*pow2(logQ3U3)*pow6(mU3)*power10(mQ3) +
         41328*msq2*pow3( logQ3U3)*pow6(mU3)*power10(mQ3) +
         716896*pow8(mU3)*power10(mQ3) + 406144*lmQ3MR*pow8(mU3)*power10(mQ3) +
         126528*dilog(1 - mQ32/mU32)* pow8(mU3)*power10(mQ3) -
         53760*lmQ3MR*dilog(1 - mQ32/mU32)*pow8(mU3)* power10(mQ3) +
         14400*dilog(1 - msq2/mU32)*pow8(mU3)*power10(mQ3) -
         157440*logQ3mq*pow8(mU3)*power10(mQ3) +
         4800*lmQ3MR*logQ3mq*pow8(mU3)*power10(mQ3) + 9600*dilog(1 -
         mQ32/mU32)*logQ3mq *pow8(mU3)*power10(mQ3) +
         390928*logQ3U3*pow8(mU3)*power10(mQ3) +
         44416*lmQ3MR*logQ3U3*pow8(mU3)*power10(mQ3) + 102816*dilog(1 -
         mQ32/mU32)*logQ3U3*pow8(mU3)*power10(mQ3) - 4320*dilog(1 - msq2/
         mU32)*logQ3U3*pow8(mU3)*power10(mQ3) -
         39840*logQ3mq*logQ3U3*pow8(mU3)*power10(mQ3) -
         7680*lmQ3MR*logQ3mq*logQ3U3*pow8(mU3)*power10(mQ3) -
         266880*pow2(lmQ3MR)*pow8(mU3)* power10(mQ3) +
         19008*logQ3U3*pow2(lmQ3MR)*pow8(mU3)*power10(mQ3) +
         9900*pow2(logQ3mq)*pow8(mU3)*power10(mQ3) -
         1890*logQ3U3*pow2(logQ3mq)*pow8(mU3)*power10(mQ3) +
         67728*pow2(logQ3U3)*pow8(mU3)*power10(mQ3) -
         62112*lmQ3MR*pow2(logQ3U3)* pow8(mU3)*power10(mQ3) +
         3360*logQ3mq*pow2(logQ3U3)*pow8( mU3)*power10(mQ3) +
         20400*pow3(logQ3U3)*pow8(mU3)*power10(mQ3) - 90*(mQ32 -
         mU32)*mU32*dilog(1 - mQ32/msq2)*(-2*mQ32 + 2*mU32 + (mQ32 +
         mU32)*logQ3U3)*(29*pow12(mQ3) + pow6(mQ3)*(-109*mU32*pow4(msq) -
         33*msq2*pow4(mU3) + 47*pow6(msq) - 13*pow6(mU3)) + pow6(msq)*pow6(mU3)
         + pow4(mQ3)*(-3*pow4(msq)*pow4(mU3) + 19*mU32*pow6(msq) +
         11*msq2*pow6( mU3)) + mQ32*(-3*pow4(mU3)*pow6(msq) +
         pow4(msq)*pow6(mU3)) + (209* msq2*mU32 - 81*pow4(msq) +
         39*pow4(mU3))*pow8(mQ3) + (5*msq2 - 119* mU32)*power10(mQ3)) -
         2880*pow4(mQ3)*pow4(msq)*power10(mU3) +
         2880*logQ3mq*pow4(mQ3)*pow4(msq)*power10(mU3) - 4320*logQ3U3*pow4(
         mQ3)*pow4(msq)*power10(mU3) + 1440*logQ3mq*logQ3U3*pow4(
         mQ3)*pow4(msq)*power10(mU3) - 450*pow2(logQ3mq)*pow4(mQ3)*pow4(
         msq)*power10(mU3) - 135*logQ3U3*pow2(logQ3mq)*pow4(mQ3)*
         pow4(msq)*power10(mU3) - 1440*pow2(logQ3U3)*pow4(mQ3)*pow4(msq)*
         power10(mU3) + 257936*msq2*pow6(mQ3)*power10(mU3) +
         249824*lmQ3MR*msq2* pow6(mQ3)*power10(mU3) + 20640*msq2*dilog(1 -
         mQ32/mU32)*pow6(mQ3)* power10(mU3) - 26880*lmQ3MR*msq2*dilog(1 -
         mQ32/mU32)*pow6(mQ3)* power10(mU3) + 14400*msq2*dilog(1 -
         msq2/mU32)*pow6(mQ3)*power10(mU3) -
         91680*msq2*logQ3mq*pow6(mQ3)*power10(mU3) + 3840*lmQ3MR*msq2*
         logQ3mq*pow6(mQ3)*power10(mU3) + 4800*msq2*dilog(1 - mQ32/mU32)*
         logQ3mq*pow6(mQ3)*power10(mU3) + 152712*msq2*logQ3U3*
         pow6(mQ3)*power10(mU3) + 151824*lmQ3MR*msq2*logQ3U3*pow6(mQ3)*
         power10(mU3) + 36528*msq2*dilog(1 - mQ32/mU32)*logQ3U3*pow6(mQ3)
         *power10(mU3) + 4320*msq2*dilog(1 - msq2/mU32)*logQ3U3*pow6(mQ3)
         *power10(mU3) - 45360*msq2*logQ3mq*logQ3U3*pow6(mQ3)* power10(mU3) -
         1920*lmQ3MR*msq2*logQ3mq*logQ3U3*pow6(mQ3) *power10(mU3) -
         141504*msq2*pow2(lmQ3MR)*pow6(mQ3)*power10(mU3) -
         25248*msq2*logQ3U3*pow2(lmQ3MR)*pow6(mQ3)*power10(mU3) + 9450*
         msq2*pow2(logQ3mq)*pow6(mQ3)*power10(mU3) +
         2835*msq2*logQ3U3*pow2(logQ3mq)*pow6(mQ3)*power10(mU3) +
         20904*msq2*pow2( logQ3U3)*pow6(mQ3)*power10(mU3) +
         2016*lmQ3MR*msq2*pow2(logQ3U3)*pow6(mQ3)*power10(mU3) -
         4320*msq2*logQ3mq*pow2(logQ3U3)*pow6(mQ3)*power10(mU3) -
         7080*msq2*pow3(logQ3U3)* pow6(mQ3)*power10(mU3) -
         450*mQ32*pow2(logQ3mq)*pow6(msq)* power10(mU3) -
         135*mQ32*logQ3U3*pow2(logQ3mq)*pow6(msq)* power10(mU3) -
         255056*pow8(mQ3)*power10(mU3) - 249824*lmQ3MR*pow8(mQ3)* power10(mU3)
         - 20640*dilog(1 - mQ32/mU32)*pow8(mQ3)*power10(mU3) +
         26880*lmQ3MR*dilog(1 - mQ32/mU32)*pow8(mQ3)*power10(mU3) -
         14400*dilog( 1 - msq2/mU32)*pow8(mQ3)*power10(mU3) +
         79200*logQ3mq*pow8(mQ3)* power10(mU3) -
         3840*lmQ3MR*logQ3mq*pow8(mQ3)*power10(mU3) - 4800*dilog(1 -
         mQ32/mU32)*logQ3mq*pow8(mQ3)*power10(mU3) -
         148392*logQ3U3*pow8(mQ3)*power10(mU3) -
         151824*lmQ3MR*logQ3U3*pow8(mQ3)*power10(mU3) - 36528*dilog(1 -
         mQ32/mU32)*logQ3U3*pow8(mQ3)*power10(mU3) - 4320*dilog(1 -
         msq2/mU32)*logQ3U3 *pow8(mQ3)*power10(mU3) +
         41040*logQ3mq*logQ3U3*pow8(mQ3) *power10(mU3) +
         1920*lmQ3MR*logQ3mq*logQ3U3*pow8(mQ3)* power10(mU3) +
         141504*pow2(lmQ3MR)*pow8(mQ3)*power10(mU3) +
         25248*logQ3U3*pow2(lmQ3MR)*pow8(mQ3)*power10(mU3) -
         8550*pow2(logQ3mq)*pow8(mQ3)*power10(mU3) - 2565*logQ3U3*pow2(logQ3mq
         )*pow8(mQ3)*power10(mU3) - 19464*pow2(logQ3U3)*pow8(mQ3)* power10(mU3)
         - 2016*lmQ3MR*pow2(logQ3U3)*pow8(mQ3)*power10(mU3) +
         4320*logQ3mq*pow2(logQ3U3)*pow8(mQ3)*power10(mU3) +
         7080*pow3(logQ3U3)*pow8(mQ3)*power10(mU3)))/(18.*(mQ32 - msq2)*
         mU32*pow4(mQ3)*pow7(mQ32 - mU32))+ (4*Xt5*(-3424*mU32*pow12(mQ3) +
         2208*lmQ3MR*mU32*pow12(mQ3) - 192*mU32* dilog(1 -
         mQ32/mU32)*pow12(mQ3) - 160*mU32*logQ3mq*pow12(mQ3) +
         3360*mU32*logQ3U3*pow12(mQ3) - 2768*lmQ3MR*mU32*logQ3U3* pow12(mQ3) +
         224*mU32*dilog(1 - mQ32/mU32)*logQ3U3*pow12(mQ3) +
         80*mU32*logQ3mq*logQ3U3*pow12(mQ3) - 256*pow14(mQ3) +
         256*lmQ3MR*pow14(mQ3) + 30*mU32*pow12(mQ3)*pow2(logQ3mq) - 15*
         mU32*logQ3U3*pow12(mQ3)*pow2(logQ3mq) - 216*mU32*pow12(
         mQ3)*pow2(logQ3U3) + 256*lmQ3MR*mU32*pow12(mQ3)*pow2(logQ3U3) -
         30*mU32*dilog(1 - msq2/mQ32)*(-2*mQ32 + 2*mU32 + (mQ32 + mU32)
         *logQ3U3)*pow2(mQ32 - msq2)*pow3(mQ32 - mU32) + 112*mU32*pow12(
         mQ3)*pow3(logQ3U3) - 1920*dilog(1 - msq2/mU32)*pow4(msq)*pow4(
         mU3)*pow6(mQ3) + 1920*logQ3mq*logQ3U3*pow4(msq)*pow4(mU3) *pow6(mQ3) -
         1080*pow2(logQ3mq)*pow4(msq)*pow4(mU3)*pow6(mQ3) +
         30*logQ3U3*pow2(logQ3mq)*pow4(msq)*pow4(mU3)*pow6(mQ3) -
         960*pow2(logQ3U3)*pow4(msq)*pow4(mU3)*pow6(mQ3) + 960*dilog(1 -
         msq2/mU32)*pow4(mQ3)*pow4(msq)*pow6(mU3) + 480*dilog(1 - msq2/mU32)*
         logQ3U3*pow4(mQ3)*pow4(msq)*pow6(mU3) -
         960*logQ3mq*logQ3U3*pow4(mQ3)*pow4(msq)*pow6(mU3) + 660*pow2(logQ3mq)*
         pow4(mQ3)*pow4(msq)*pow6(mU3) + 240*logQ3U3*pow2(logQ3mq)
         *pow4(mQ3)*pow4(msq)*pow6(mU3) + 480*pow2(logQ3U3)*pow4(mQ3)*
         pow4(msq)*pow6(mU3) - 480*logQ3mq*pow2(logQ3U3)*pow4(mQ3)
         *pow4(msq)*pow6(mU3) + 240*pow3(logQ3U3)*pow4(mQ3)*pow4(msq)*
         pow6(mU3) + 2880*msq2*pow6(mQ3)*pow6(mU3) + 3840*msq2*dilog(1 - msq2/
         mU32)*pow6(mQ3)*pow6(mU3) - 2880*msq2*logQ3mq*pow6(mQ3)*pow6( mU3) +
         2400*msq2*logQ3U3*pow6(mQ3)*pow6(mU3) -
         4320*msq2*logQ3mq*logQ3U3*pow6(mQ3)*pow6(mU3) +
         1560*msq2*pow2(logQ3mq)*pow6(mQ3)*pow6(mU3) +
         1920*msq2*pow2(logQ3U3)*pow6( mQ3)*pow6(mU3) + 960*mU32*dilog(1 -
         msq2/mU32)*pow4(msq)*pow8(mQ3) - 480*mU32*dilog(1 -
         msq2/mU32)*logQ3U3*pow4(msq)*pow8(mQ3) - 960*
         mU32*logQ3mq*logQ3U3*pow4(msq)*pow8(mQ3) + 510*mU32*pow2(
         logQ3mq)*pow4(msq)*pow8(mQ3) -
         255*mU32*logQ3U3*pow2(logQ3mq)*pow4(msq)*pow8(mQ3) +
         480*mU32*pow2(logQ3U3)*pow4( msq)*pow8(mQ3) +
         480*mU32*logQ3mq*pow2(logQ3U3)*pow4(msq) *pow8(mQ3) -
         240*mU32*pow3(logQ3U3)*pow4(msq)*pow8(mQ3) - 2880*
         msq2*pow4(mU3)*pow8(mQ3) - 1920*msq2*dilog(1 - msq2/mU32)*pow4(mU3)*
         pow8(mQ3) + 2880*msq2*logQ3mq*pow4(mU3)*pow8(mQ3) - 480*msq2*
         logQ3U3*pow4(mU3)*pow8(mQ3) + 960*msq2*dilog(1 -
         msq2/mU32)*logQ3U3*pow4(mU3)*pow8(mQ3) +
         1440*msq2*logQ3mq*logQ3U3*pow4(mU3)*pow8(mQ3) -
         720*msq2*pow2(logQ3mq)*pow4(mU3)* pow8(mQ3) +
         420*msq2*logQ3U3*pow2(logQ3mq)*pow4(mU3)* pow8(mQ3) -
         480*msq2*pow2(logQ3U3)*pow4(mU3)*pow8(mQ3) - 960*
         msq2*logQ3mq*pow2(logQ3U3)*pow4(mU3)*pow8(mQ3) + 480*
         msq2*pow3(logQ3U3)*pow4(mU3)*pow8(mQ3) - 24128*pow6(mU3)*pow8( mQ3) +
         18560*lmQ3MR*pow6(mU3)*pow8(mQ3) - 1152*dilog(1 - mQ32/mU32)*
         pow6(mU3)*pow8(mQ3) + 960*dilog(1 - msq2/mU32)*pow6(mU3)*pow8(mQ3) -
         960*logQ3mq*pow6(mU3)*pow8(mQ3) - 12160*logQ3U3*pow6(mU3) *pow8(mQ3) +
         4256*lmQ3MR*logQ3U3*pow6(mU3)*pow8(mQ3) + 768* dilog(1 -
         mQ32/mU32)*logQ3U3*pow6(mU3)*pow8(mQ3) - 480*dilog(1 -
         msq2/mU32)*logQ3U3*pow6(mU3)*pow8(mQ3) -
         480*logQ3mq*logQ3U3*pow6(mU3)*pow8(mQ3) + 660*pow2(logQ3mq)*pow6(mU3)*
         pow8(mQ3) - 240*logQ3U3*pow2(logQ3mq)*pow6(mU3)*pow8(mQ3) +
         704*pow2(logQ3U3)*pow6(mU3)*pow8(mQ3) -
         2528*lmQ3MR*pow2(logQ3U3)*pow6(mU3)*pow8(mQ3) +
         80*logQ3mq*pow2(logQ3U3 )*pow6(mU3)*pow8(mQ3) +
         1352*pow3(logQ3U3)*pow6(mU3)*pow8(mQ3) - 30*(mQ32 - mU32)*mU32*dilog(1
         - mQ32/msq2)*(-2*mQ32 + 2*mU32 + (mQ32 +
         mU32)*logQ3U3)*(-2*mQ32*msq2*mU32*(msq2 + mU32) + pow4(msq)* pow4(mU3)
         + pow4(mQ3)*(-28*msq2*mU32 + 17*pow4(msq) + 17*pow4(mU3)) - 2*(msq2 +
         mU32)*pow6(mQ3) + pow8(mQ3)) - 960*msq2*pow4(mQ3)*pow8(mU3) -
         1920*msq2*dilog(1 - msq2/mU32)*pow4(mQ3)*pow8(mU3) +
         960*msq2*logQ3mq*pow4(mQ3)*pow8(mU3) -
         1440*msq2*logQ3U3*pow4(mQ3)*pow8( mU3) - 960*msq2*dilog(1 -
         msq2/mU32)*logQ3U3*pow4(mQ3)*pow8(mU3) +
         2400*msq2*logQ3mq*logQ3U3*pow4(mQ3)*pow8(mU3) - 720*
         msq2*pow2(logQ3mq)*pow4(mQ3)*pow8(mU3) - 420*msq2*logQ3U3
         *pow2(logQ3mq)*pow4(mQ3)*pow8(mU3) -
         1440*msq2*pow2(logQ3U3)*pow4(mQ3)*pow8(mU3) +
         960*msq2*logQ3mq*pow2(logQ3U3)*pow4(mQ3)*pow8(mU3) -
         480*msq2*pow3(logQ3U3)*pow4(mQ3)* pow8(mU3) -
         120*mQ32*pow2(logQ3mq)*pow4(msq)*pow8(mU3) - 30*
         mQ32*logQ3U3*pow2(logQ3mq)*pow4(msq)*pow8(mU3) + 16000*
         pow6(mQ3)*pow8(mU3) - 12864*lmQ3MR*pow6(mQ3)*pow8(mU3) + 768*dilog(1 -
         mQ32/mU32)*pow6(mQ3)*pow8(mU3) - 1920*dilog(1 - msq2/mU32)*pow6(mQ3)*
         pow8(mU3) + 640*logQ3mq*pow6(mQ3)*pow8(mU3) +
         17600*logQ3U3*pow6(mQ3)*pow8(mU3) -
         9920*lmQ3MR*logQ3U3*pow6(mQ3)*pow8( mU3) - 320*dilog(1 -
         mQ32/mU32)*logQ3U3*pow6(mQ3)*pow8(mU3) +
         640*logQ3mq*logQ3U3*pow6(mQ3)*pow8(mU3) -
         1080*pow2(logQ3mq)*pow6(mQ3)*pow8(mU3) - 30*logQ3U3*pow2(logQ3mq
         )*pow6(mQ3)*pow8(mU3) + 5872*pow2(logQ3U3)*pow6(mQ3)*pow8(mU3) -
         448*lmQ3MR*pow2(logQ3U3)*pow6(mQ3)*pow8(mU3) + 80*logQ3mq
         *pow2(logQ3U3)*pow6(mQ3)*pow8(mU3) + 232*pow3(logQ3U3)*
         pow6(mQ3)*pow8(mU3) + 960*msq2*mU32*power10(mQ3) -
         960*msq2*mU32*logQ3mq*power10(mQ3) -
         480*msq2*mU32*logQ3U3*power10(mQ3) +
         480*msq2*mU32*logQ3mq*logQ3U3*power10(mQ3) - 60*msq2*
         mU32*pow2(logQ3mq)*power10(mQ3) + 30*msq2*mU32*logQ3U3*
         pow2(logQ3mq)*power10(mQ3) + 15744*pow4(mU3)*power10(mQ3) -
         11456*lmQ3MR*pow4(mU3)*power10(mQ3) + 768*dilog(1 - mQ32/mU32)*pow4(
         mU3)*power10(mQ3) + 640*logQ3mq*pow4(mU3)*power10(mQ3) - 2240*
         logQ3U3*pow4(mU3)*power10(mQ3) + 4160*lmQ3MR*logQ3U3*
         pow4(mU3)*power10(mQ3) - 704*dilog(1 - mQ32/mU32)*logQ3U3*pow4(
         mU3)*power10(mQ3) - 120*pow2(logQ3mq)*pow4(mU3)*power10(mQ3) +
         30*logQ3U3*pow2(logQ3mq)*pow4(mU3)*power10(mQ3) - 2480*
         pow2(logQ3U3)*pow4(mU3)*power10(mQ3) +
         1472*lmQ3MR*pow2(logQ3U3)*pow4(mU3)*power10(mQ3) -
         80*logQ3mq*pow2(logQ3U3)*pow4(mU3)*power10(mQ3) -
         616*pow3(logQ3U3)*pow4(mU3)* power10(mQ3) -
         60*mQ32*msq2*pow2(logQ3mq)*power10(mU3) - 30*
         mQ32*msq2*logQ3U3*pow2(logQ3mq)*power10(mU3) - 3936*pow4(
         mQ3)*power10(mU3) + 3296*lmQ3MR*pow4(mQ3)*power10(mU3) - 192*dilog(1 -
         mQ32/mU32)*pow4(mQ3)*power10(mU3) + 960*dilog(1 -
         msq2/mU32)*pow4(mQ3)* power10(mU3) -
         160*logQ3mq*pow4(mQ3)*power10(mU3) -
         6560*logQ3U3*pow4(mQ3)*power10(mU3) + 4272*lmQ3MR*logQ3U3*pow4(
         mQ3)*power10(mU3) + 32*dilog(1 - mQ32/mU32)*logQ3U3*pow4(mQ3)*
         power10(mU3) + 480*dilog(1 - msq2/mU32)*logQ3U3*pow4(mQ3)*
         power10(mU3) - 240*logQ3mq*logQ3U3*pow4(mQ3)*power10(mU3) +
         510*pow2(logQ3mq)*pow4(mQ3)*power10(mU3) + 255*logQ3U3*
         pow2(logQ3mq)*pow4(mQ3)*power10(mU3) - 3880*pow2(logQ3U3)
         *pow4(mQ3)*power10(mU3) + 1248*lmQ3MR*pow2(logQ3U3)*pow4(mQ3)*
         power10(mU3) - 80*logQ3mq*pow2(logQ3U3)*pow4(mQ3)* power10(mU3) -
         824*pow3(logQ3U3)*pow4(mQ3)*power10(mU3) + 30*
         pow2(logQ3mq)*pow4(msq)*power10(mU3) + 15*logQ3U3*pow2(
         logQ3mq)*pow4(msq)*power10(mU3)))/(3.*mU32*pow3(mQ3)*pow7(mQ32 -
         mU32))+ (1280*mQ32*Xt6*(-2*mQ32 + 2*mU32 + (mQ32 +
         mU32)*logQ3U3)*pow2(- mQ32 + mU32 + mU32*logQ3U3))/(3.*pow7(mQ32 -
         mU32));
      }

      case Limits::MU3_EQ_M3: {
         const double logQ3U3 = std::log(mQ32/mU32);
         const double logU3Q3 = -logQ3U3;
         const double logQ3mq = std::log(mQ32/msq2);
         const double logU3mq = std::log(mU32/msq2);

         return
         (Xt4*(84*msq2*mU32*dilog(1 - mU32/mQ32)*pow14(mQ3) -
         78*msq2*mU32*dilog(1 - mU32/mQ32)*logQ3U3*pow14(mQ3) +
         244640*mQ32*msq2*pow14(mU3) - 13120*lmQ3MR*mQ32*msq2*pow14(mU3) +
         1800*mQ32*msq2*dilog(1 - msq2/mU32) *pow14(mU3) +
         80556*mQ32*msq2*dilog(1 - mU32/mQ32)*pow14(mU3) - 10752*
         lmQ3MR*mQ32*msq2*dilog(1 - mU32/mQ32)*pow14(mU3) + 1800*mQ32*msq2*
         dilog(1 - mU32/msq2)*pow14(mU3) - 10560*mQ32*msq2*logQ3mq*pow14( mU3)
         - 7680*lmQ3MR*mQ32*msq2*logQ3mq*pow14(mU3) + 1920*mQ32* msq2*dilog(1 -
         mU32/mQ32)*logQ3mq*pow14(mU3) + 114224*mQ32*msq2* logQ3U3*pow14(mU3) +
         70240*lmQ3MR*mQ32*msq2*logQ3U3* pow14(mU3) + 900*mQ32*msq2*dilog(1 -
         msq2/mU32)*logQ3U3*pow14( mU3) + 24534*mQ32*msq2*dilog(1 -
         mU32/mQ32)*logQ3U3*pow14(mU3) + 900*mQ32*msq2*dilog(1 -
         mU32/msq2)*logQ3U3*pow14(mU3) - 7080*
         mQ32*msq2*logQ3mq*logQ3U3*pow14(mU3) + 3840*lmQ3MR*mQ32*
         msq2*logQ3mq*logQ3U3*pow14(mU3) - 384*mQ32*msq2*logU3Q3*pow14(mU3) -
         192*mQ32*msq2*logQ3U3*logU3Q3*pow14( mU3) + 36*msq2*dilog(1 -
         mU32/mQ32)*pow16(mQ3) - 36*mU32*dilog(1 - mU32/mQ32)*pow16(mQ3) -
         18*msq2*dilog(1 - mU32/mQ32)*logQ3U3* pow16(mQ3) + 18*mU32*dilog(1 -
         mU32/mQ32)*logQ3U3*pow16(mQ3) - 261920*mQ32*pow16(mU3) +
         13120*lmQ3MR*mQ32*pow16(mU3) + 384*msq2*pow16( mU3) +
         10368*lmQ3MR*msq2*pow16(mU3) - 1080*mQ32*dilog(1 - msq2/mU32)*
         pow16(mU3) - 80556*mQ32*dilog(1 - mU32/mQ32)*pow16(mU3) +
         10752*lmQ3MR* mQ32*dilog(1 - mU32/mQ32)*pow16(mU3) +
         10440*mQ32*dilog(1 - mU32/msq2)* pow16(mU3) +
         27840*mQ32*logQ3mq*pow16(mU3) + 7680*lmQ3MR*mQ32* logQ3mq*pow16(mU3) -
         1920*msq2*logQ3mq*pow16(mU3) + 1920* lmQ3MR*msq2*logQ3mq*pow16(mU3) -
         1920*mQ32*dilog(1 - mU32/mQ32)* logQ3mq*pow16(mU3) -
         140144*mQ32*logQ3U3*pow16(mU3) - 70240*lmQ3MR*mQ32*logQ3U3*pow16(mU3)
         - 8256*msq2*logQ3U3* pow16(mU3) + 19392*lmQ3MR*msq2*logQ3U3*pow16(mU3)
         - 540*mQ32* dilog(1 - msq2/mU32)*logQ3U3*pow16(mU3) -
         24534*mQ32*dilog(1 - mU32/mQ32)*logQ3U3*pow16(mU3) + 5220*mQ32*dilog(1
         - mU32/msq2)* logQ3U3*pow16(mU3) + 15000*mQ32*logQ3mq*logQ3U3*
         pow16(mU3) - 3840*lmQ3MR*mQ32*logQ3mq*logQ3U3*pow16(mU3) -
         1920*msq2*logQ3mq*logQ3U3*pow16(mU3) + 384*mQ32*logU3Q3*pow16(mU3) +
         192*mQ32*logQ3U3*logU3Q3*pow16( mU3) - 3840*mQ32*logU3mq*pow16(mU3) -
         1920*mQ32*logQ3U3* logU3mq*pow16(mU3) - 384*pow18(mU3) -
         10368*lmQ3MR*pow18(mU3) + 1920*logQ3mq*pow18(mU3) -
         1920*lmQ3MR*logQ3mq*pow18(mU3) + 8256*logQ3U3*pow18(mU3) -
         19392*lmQ3MR*logQ3U3*pow18( mU3) + 1920*logQ3mq*logQ3U3*pow18(mU3) -
         4992*mQ32*msq2* pow14(mU3)*pow2(lmQ3MR) -
         45504*mQ32*msq2*logQ3U3*pow14(mU3)* pow2(lmQ3MR) +
         4992*mQ32*pow16(mU3)*pow2(lmQ3MR) - 10752*msq2*pow16(
         mU3)*pow2(lmQ3MR) + 45504*mQ32*logQ3U3*pow16(mU3)*pow2(lmQ3MR) +
         10752*pow18(mU3)*pow2(lmQ3MR) +
         6660*mQ32*msq2*pow14(mU3)*pow2(logQ3mq) +
         3330*mQ32*msq2*logQ3U3*pow14(mU3)*pow2(logQ3mq) -
         540*mQ32*pow16(mU3)*pow2(logQ3mq) -
         270*mQ32*logQ3U3*pow16(mU3)*pow2(logQ3mq) +
         42*msq2*mU32*pow14(mQ3)*pow2( logQ3U3) -
         15094*mQ32*msq2*pow14(mU3)*pow2(logQ3U3) +
         82176*lmQ3MR*mQ32*msq2*pow14(mU3)*pow2(logQ3U3) - 2820*mQ32*
         msq2*logQ3mq*pow14(mU3)*pow2(logQ3U3) + 18*msq2*pow16(
         mQ3)*pow2(logQ3U3) - 18*mU32*pow16(mQ3)*pow2(logQ3U3) +
         6814*mQ32*pow16(mU3)*pow2(logQ3U3) - 82176*lmQ3MR*mQ32*pow16(
         mU3)*pow2(logQ3U3) - 8640*msq2*pow16(mU3)*pow2(logQ3U3) +
         2460*mQ32*logQ3mq*pow16(mU3)*pow2(logQ3U3) + 8640*pow18(
         mU3)*pow2(logQ3U3) - 39*msq2*mU32*pow14(mQ3)*pow3(logQ3U3 ) -
         30099*mQ32*msq2*pow14(mU3)*pow3(logQ3U3) - 9*msq2*pow16(mQ3)
         *pow3(logQ3U3) + 9*mU32*pow16(mQ3)*pow3(logQ3U3) + 30279*
         mQ32*pow16(mU3)*pow3(logQ3U3) - 1072288*msq2*pow12(mU3)*pow4( mQ3) -
         170176*lmQ3MR*msq2*pow12(mU3)*pow4(mQ3) - 9000*msq2*dilog(1 -
         msq2/mU32)*pow12(mU3)*pow4(mQ3) - 305508*msq2*dilog(1 - mU32/mQ32)*
         pow12(mU3)*pow4(mQ3) + 53760*lmQ3MR*msq2*dilog(1 - mU32/mQ32)*pow12(
         mU3)*pow4(mQ3) + 71640*msq2*dilog(1 - mU32/msq2)*pow12(mU3)*pow4(mQ3)
         + 120000*msq2*logQ3mq*pow12(mU3)*pow4(mQ3) +
         9600*lmQ3MR*msq2*logQ3mq*pow12(mU3)*pow4(mQ3) - 9600*msq2*dilog(1 -
         mU32/mQ32)*logQ3mq*pow12(mU3)*pow4(mQ3) +
         24880*msq2*logQ3U3*pow12(mU3)* pow4(mQ3) -
         621600*lmQ3MR*msq2*logQ3U3*pow12(mU3)*pow4(mQ3) - 2700*msq2*dilog(1 -
         msq2/mU32)*logQ3U3*pow12(mU3)*pow4(mQ3) - 13062*msq2*dilog(1 -
         mU32/mQ32)*logQ3U3*pow12(mU3)*pow4(mQ3) + 37620*msq2*dilog(1 -
         mU32/msq2)*logQ3U3*pow12(mU3)*pow4(mQ3) +
         38280*msq2*logQ3mq*logQ3U3*pow12(mU3)*pow4(mQ3) - 15360*
         lmQ3MR*msq2*logQ3mq*logQ3U3*pow12(mU3)*pow4(mQ3) + 768*
         msq2*logU3Q3*pow12(mU3)*pow4(mQ3) + 1118368*pow14(mU3)*pow4(mQ3) +
         170176*lmQ3MR*pow14(mU3)*pow4(mQ3) + 5400*dilog(1 - msq2/mU32)*pow14(
         mU3)*pow4(mQ3) + 305508*dilog(1 - mU32/mQ32)*pow14(mU3)*pow4(mQ3) -
         53760*lmQ3MR*dilog(1 - mU32/mQ32)*pow14(mU3)*pow4(mQ3) - 63720*dilog(1
         - mU32/msq2)*pow14(mU3)*pow4(mQ3) - 166080*logQ3mq*pow14(mU3)*
         pow4(mQ3) - 9600*lmQ3MR*logQ3mq*pow14(mU3)*pow4(mQ3) + 9600* dilog(1 -
         mU32/mQ32)*logQ3mq*pow14(mU3)*pow4(mQ3) +
         9680*logQ3U3*pow14(mU3)*pow4(mQ3) + 621600*lmQ3MR*logQ3U3*pow14(
         mU3)*pow4(mQ3) + 1620*dilog(1 - msq2/mU32)*logQ3U3*pow14(mU3)*
         pow4(mQ3) + 13062*dilog(1 - mU32/mQ32)*logQ3U3*pow14(mU3)*pow4( mQ3) -
         21420*dilog(1 - mU32/msq2)*logQ3U3*pow14(mU3)*pow4(mQ3) -
         40440*logQ3mq*logQ3U3*pow14(mU3)*pow4(mQ3) + 15360*
         lmQ3MR*logQ3mq*logQ3U3*pow14(mU3)*pow4(mQ3) -
         768*logU3Q3*pow14(mU3)*pow4(mQ3) + 19200*logU3mq*pow14(mU3)*pow4( mQ3)
         + 5760*logQ3U3*logU3mq*pow14(mU3)*pow4(mQ3) +
         186240*msq2*pow12(mU3)*pow2(lmQ3MR)*pow4(mQ3) +
         158016*msq2*logQ3U3*pow12(mU3)*pow2(lmQ3MR)*pow4(mQ3) -
         186240*pow14(mU3)*pow2( lmQ3MR)*pow4(mQ3) -
         158016*logQ3U3*pow14(mU3)*pow2(lmQ3MR)*pow4( mQ3) +
         7020*msq2*pow12(mU3)*pow2(logQ3mq)*pow4(mQ3) + 10170*
         msq2*logQ3U3*pow12(mU3)*pow2(logQ3mq)*pow4(mQ3) - 3060*
         pow14(mU3)*pow2(logQ3mq)*pow4(mQ3) - 2070*logQ3U3*pow14(
         mU3)*pow2(logQ3mq)*pow4(mQ3) +
         327498*msq2*pow12(mU3)*pow2(logQ3U3)*pow4(mQ3) -
         264192*lmQ3MR*msq2*pow12(mU3)*pow2(logQ3U3)*pow4(mQ3) -
         1140*msq2*logQ3mq*pow12(mU3)*pow2(logQ3U3)*pow4(mQ3) -
         332178*pow14(mU3)*pow2(logQ3U3)*pow4(mQ3) +
         264192*lmQ3MR*pow14(mU3)*pow2(logQ3U3)*pow4(mQ3) +
         2220*logQ3mq*pow14(mU3)*pow2(logQ3U3)*pow4(mQ3) + 145527*msq2*
         pow12(mU3)*pow3(logQ3U3)*pow4(mQ3) -
         146067*pow14(mU3)*pow3(logQ3U3)*pow4(mQ3) + 360*mU32*dilog(1 -
         msq2/mU32)*pow12(mQ3)*pow4( msq) + 360*mU32*dilog(1 -
         mU32/msq2)*pow12(mQ3)*pow4(msq) - 180*mU32* dilog(1 -
         msq2/mU32)*logQ3U3*pow12(mQ3)*pow4(msq) - 180*mU32* dilog(1 -
         mU32/msq2)*logQ3U3*pow12(mQ3)*pow4(msq) - 360*mU32*
         logQ3mq*logQ3U3*pow12(mQ3)*pow4(msq) + 17280*mQ32*pow12(
         mU3)*pow4(msq) - 360*mQ32*dilog(1 - msq2/mU32)*pow12(mU3)*pow4(msq) -
         29160*mQ32*dilog(1 - mU32/msq2)*pow12(mU3)*pow4(msq) -
         17280*mQ32*logQ3mq*pow12(mU3)*pow4(msq) +
         25920*mQ32*logQ3U3*pow12(mU3)* pow4(msq) - 180*mQ32*dilog(1 -
         msq2/mU32)*logQ3U3*pow12(mU3)* pow4(msq) - 14580*mQ32*dilog(1 -
         mU32/msq2)*logQ3U3*pow12(mU3)* pow4(msq) -
         8280*mQ32*logQ3mq*logQ3U3*pow12(mU3)*pow4( msq) +
         180*mU32*pow12(mQ3)*pow2(logQ3mq)*pow4(msq) - 90*mU32*
         logQ3U3*pow12(mQ3)*pow2(logQ3mq)*pow4(msq) - 14580*mQ32*
         pow12(mU3)*pow2(logQ3mq)*pow4(msq) - 7290*mQ32*logQ3U3*
         pow12(mU3)*pow2(logQ3mq)*pow4(msq) + 180*mU32*pow12(mQ3)*pow2(
         logQ3U3)*pow4(msq) + 180*mU32*logQ3mq*pow12(mQ3)*pow2(
         logQ3U3)*pow4(msq) + 8460*mQ32*pow12(mU3)*pow2(logQ3U3)* pow4(msq) +
         180*mQ32*logQ3mq*pow12(mU3)*pow2(logQ3U3)* pow4(msq) -
         90*mU32*pow12(mQ3)*pow3(logQ3U3)*pow4(msq) - 90*
         mQ32*pow12(mU3)*pow3(logQ3U3)*pow4(msq) - 59936*msq2*pow12(mQ3)*
         pow4(mU3) - 117440*lmQ3MR*msq2*pow12(mQ3)*pow4(mU3) -
         1800*msq2*dilog(1 - msq2/mU32)*pow12(mQ3)*pow4(mU3) +
         6804*msq2*dilog(1 - mU32/mQ32)* pow12(mQ3)*pow4(mU3) +
         10752*lmQ3MR*msq2*dilog(1 - mU32/mQ32)*pow12( mQ3)*pow4(mU3) +
         3960*msq2*dilog(1 - mU32/msq2)*pow12(mQ3)*pow4(mU3) +
         35520*msq2*logQ3mq*pow12(mQ3)*pow4(mU3) -
         1920*lmQ3MR*msq2*logQ3mq*pow12(mQ3)*pow4(mU3) - 1920*msq2*dilog(1 -
         mU32/mQ32)*logQ3mq*pow12(mQ3)*pow4(mU3) +
         80112*msq2*logQ3U3*pow12(mQ3)* pow4(mU3) +
         12256*lmQ3MR*msq2*logQ3U3*pow12(mQ3)*pow4(mU3) + 900*msq2*dilog(1 -
         msq2/mU32)*logQ3U3*pow12(mQ3)*pow4(mU3) - 7362*msq2*dilog(1 -
         mU32/mQ32)*logQ3U3*pow12(mQ3)*pow4(mU3) - 1980*msq2*dilog(1 -
         mU32/msq2)*logQ3U3*pow12(mQ3)*pow4(mU3) -
         16920*msq2*logQ3mq*logQ3U3*pow12(mQ3)*pow4(mU3) - 768*
         msq2*logU3Q3*pow12(mQ3)*pow4(mU3) +
         384*msq2*logQ3U3*logU3Q3*pow12(mQ3)*pow4(mU3) - 84*dilog(1 -
         mU32/mQ32)*pow14(mQ3)* pow4(mU3) + 78*dilog(1 -
         mU32/mQ32)*logQ3U3*pow14(mQ3)*pow4(mU3) +
         58752*msq2*pow12(mQ3)*pow2(lmQ3MR)*pow4(mU3) -
         24000*msq2*logQ3U3*pow12(mQ3)*pow2(lmQ3MR)*pow4(mU3) -
         3780*msq2*pow12(mQ3)*pow2( logQ3mq)*pow4(mU3) +
         1890*msq2*logQ3U3*pow12(mQ3)*pow2( logQ3mq)*pow4(mU3) -
         19306*msq2*pow12(mQ3)*pow2(logQ3U3)* pow4(mU3) +
         13248*lmQ3MR*msq2*pow12(mQ3)*pow2(logQ3U3)*pow4(mU3) +
         1020*msq2*logQ3mq*pow12(mQ3)*pow2(logQ3U3)*pow4(mU3) -
         42*pow14(mQ3)*pow2(logQ3U3)*pow4(mU3) - 63*msq2*pow12(mQ3)*pow3(
         logQ3U3)*pow4(mU3) + 39*pow14(mQ3)*pow3(logQ3U3)*pow4( mU3) +
         2880*mQ32*(mQ32 - msq2)*(mQ32 - mU32)*(msq2 - mU32)*dilog(1 -
         msq2/mQ32)*(-2*mQ32 + 2*mU32 + (mQ32 + mU32)*logQ3U3)*(mQ32*( msq2 -
         3*mU32) + 3*msq2*mU32 + pow4(mQ3) - 2*pow4(mU3))*pow4(mU3) -
         1839680*pow12(mU3)*pow6(mQ3) - 603008*lmQ3MR*pow12(mU3)*pow6(mQ3) -
         10800*dilog(1 - msq2/mU32)*pow12(mU3)*pow6(mQ3) - 425916*dilog(1 -
         mU32/mQ32)*pow12(mU3)*pow6(mQ3) + 107520*lmQ3MR*dilog(1 - mU32/mQ32)*
         pow12(mU3)*pow6(mQ3) + 110160*dilog(1 -
         mU32/msq2)*pow12(mU3)*pow6(mQ3) + 347520*logQ3mq*pow12(mU3)*pow6(mQ3)
         - 19200*dilog(1 - mU32/ mQ32)*logQ3mq*pow12(mU3)*pow6(mQ3) +
         823968*logQ3U3* pow12(mU3)*pow6(mQ3) -
         1234112*lmQ3MR*logQ3U3*pow12(mU3)*pow6( mQ3) - 1080*dilog(1 -
         msq2/mU32)*logQ3U3*pow12(mU3)*pow6(mQ3) + 91098*dilog(1 -
         mU32/mQ32)*logQ3U3*pow12(mU3)*pow6(mQ3) + 1800* dilog(1 -
         mU32/msq2)*logQ3U3*pow12(mU3)*pow6(mQ3) +
         17520*logQ3mq*logQ3U3*pow12(mU3)*pow6(mQ3) -
         23040*lmQ3MR*logQ3mq*logQ3U3*pow12(mU3)*pow6(mQ3) - 768*logU3Q3*pow12(
         mU3)*pow6(mQ3) - 768*logQ3U3*logU3Q3*pow12(mU3)*pow6(mQ3) -
         38400*logU3mq*pow12(mU3)*pow6(mQ3) -
         3840*logQ3U3*logU3mq*pow12(mU3)*pow6(mQ3) +
         480000*pow12(mU3)*pow2(lmQ3MR)*pow6( mQ3) +
         177024*logQ3U3*pow12(mU3)*pow2(lmQ3MR)*pow6(mQ3) - 2520*
         pow12(mU3)*pow2(logQ3mq)*pow6(mQ3) - 4860*logQ3U3*pow12(
         mU3)*pow2(logQ3mq)*pow6(mQ3) +
         658922*pow12(mU3)*pow2(logQ3U3)*pow6(mQ3) -
         206016*lmQ3MR*pow12(mU3)*pow2(logQ3U3)*pow6( mQ3) +
         10680*logQ3mq*pow12(mU3)*pow2(logQ3U3)*pow6(mQ3) +
         125745*pow12(mU3)*pow3(logQ3U3)*pow6(mQ3) + 360*dilog(1 - msq2/
         mU32)*pow12(mQ3)*pow6(msq) + 360*dilog(1 - mU32/msq2)*pow12(mQ3)*pow6(
         msq) - 180*dilog(1 - msq2/mU32)*logQ3U3*pow12(mQ3)*pow6(msq) -
         180*dilog(1 - mU32/msq2)*logQ3U3*pow12(mQ3)*pow6(msq) -
         360*logQ3mq*logQ3U3*pow12(mQ3)*pow6(msq) + 180*pow12(mQ3)*pow2(
         logQ3mq)*pow6(msq) - 90*logQ3U3*pow12(mQ3)*pow2(logQ3mq)*pow6(msq) +
         180*pow12(mQ3)*pow2(logQ3U3)*pow6(msq) + 180*
         logQ3mq*pow12(mQ3)*pow2(logQ3U3)*pow6(msq) - 90*pow12(
         mQ3)*pow3(logQ3U3)*pow6(msq) + 6*mQ32*(-msq2 + mU32)*dilog(1 -
         mQ32/mU32)*(-2*mQ32 + 2*mU32 + (mQ32 + mU32)*logQ3U3)*pow3(mQ32 -
         mU32)*(19*mU32*pow4(mQ3) - 31*mQ32*pow4(mU3) + 3*pow6(mQ3) - 119*
         pow6(mU3)) + 59936*pow12(mQ3)*pow6(mU3) + 117440*lmQ3MR*pow12(mQ3)*
         pow6(mU3) + 1080*dilog(1 - msq2/mU32)*pow12(mQ3)*pow6(mU3) - 6804*
         dilog(1 - mU32/mQ32)*pow12(mQ3)*pow6(mU3) - 10752*lmQ3MR*dilog(1 -
         mU32/mQ32)*pow12(mQ3)*pow6(mU3) - 4680*dilog(1 -
         mU32/msq2)*pow12(mQ3)* pow6(mU3) - 35520*logQ3mq*pow12(mQ3)*pow6(mU3)
         + 1920*lmQ3MR* logQ3mq*pow12(mQ3)*pow6(mU3) + 1920*dilog(1 -
         mU32/mQ32)*logQ3mq*pow12(mQ3)*pow6(mU3) -
         80112*logQ3U3*pow12(mQ3)*pow6( mU3) -
         12256*lmQ3MR*logQ3U3*pow12(mQ3)*pow6(mU3) - 540*dilog(1 -
         msq2/mU32)*logQ3U3*pow12(mQ3)*pow6(mU3) + 7362*dilog(1 - mU32/
         mQ32)*logQ3U3*pow12(mQ3)*pow6(mU3) + 2340*dilog(1 - mU32/msq2)*
         logQ3U3*pow12(mQ3)*pow6(mU3) +
         17640*logQ3mq*logQ3U3*pow12(mQ3)*pow6(mU3) +
         768*logU3Q3*pow12(mQ3)*pow6(mU3) -
         384*logQ3U3*logU3Q3*pow12(mQ3)*pow6(mU3) +
         3840*logU3mq*pow12(mQ3)*pow6(mU3) - 1920*logQ3U3*logU3mq*pow12(
         mQ3)*pow6(mU3) - 58752*pow12(mQ3)*pow2(lmQ3MR)*pow6(mU3) +
         24000*logQ3U3*pow12(mQ3)*pow2(lmQ3MR)*pow6(mU3) +
         3420*pow12(mQ3)*pow2( logQ3mq)*pow6(mU3) -
         1710*logQ3U3*pow12(mQ3)*pow2(logQ3mq)*pow6(mU3) +
         18946*pow12(mQ3)*pow2(logQ3U3)*pow6(mU3) -
         13248*lmQ3MR*pow12(mQ3)*pow2(logQ3U3)*pow6(mU3) -
         1380*logQ3mq*pow12(mQ3)*pow2(logQ3U3)*pow6(mU3) + 243*pow12(mQ3)*
         pow3(logQ3U3)*pow6(mU3) - 3600*dilog(1 - msq2/mU32)*pow6(mQ3)*
         pow6(msq)*pow6(mU3) + 2160*dilog(1 - mU32/msq2)*pow6(mQ3)*pow6(msq)*
         pow6(mU3) - 360*dilog(1 - msq2/mU32)*logQ3U3*pow6(mQ3)*pow6(msq)
         *pow6(mU3) - 9000*dilog(1 - mU32/msq2)*logQ3U3*pow6(mQ3)*pow6(
         msq)*pow6(mU3) + 3600*logQ3mq*logQ3U3*pow6(mQ3)*pow6(msq) *pow6(mU3) +
         1080*pow2(logQ3mq)*pow6(mQ3)*pow6(msq)*pow6(mU3) -
         4500*logQ3U3*pow2(logQ3mq)*pow6(mQ3)*pow6(msq)*pow6(mU3) -
         1800*pow2(logQ3U3)*pow6(mQ3)*pow6(msq)*pow6(mU3) +
         360*logQ3mq*pow2(logQ3U3)*pow6(mQ3)*pow6(msq)*pow6(mU3) - 180*
         pow3(logQ3U3)*pow6(mQ3)*pow6(msq)*pow6(mU3) + 3600*dilog(1 -
         msq2/mU32)*pow4(mU3)*pow6(msq)*pow8(mQ3) + 9360*dilog(1 - mU32/msq2)*
         pow4(mU3)*pow6(msq)*pow8(mQ3) - 360*dilog(1 - msq2/mU32)*logQ3U3
         *pow4(mU3)*pow6(msq)*pow8(mQ3) - 3240*dilog(1 -
         mU32/msq2)*logQ3U3*pow4(mU3)*pow6(msq)*pow8(mQ3) -
         3600*logQ3mq*logQ3U3*pow4(mU3)*pow6(msq)*pow8(mQ3) +
         4680*pow2(logQ3mq)*pow4( mU3)*pow6(msq)*pow8(mQ3) -
         1620*logQ3U3*pow2(logQ3mq)* pow4(mU3)*pow6(msq)*pow8(mQ3) +
         1800*pow2(logQ3U3)*pow4(mU3)* pow6(msq)*pow8(mQ3) +
         360*logQ3mq*pow2(logQ3U3)*pow4(mU3) *pow6(msq)*pow8(mQ3) -
         180*pow3(logQ3U3)*pow4(mU3)*pow6(msq)* pow8(mQ3) + 3600*dilog(1 -
         msq2/mU32)*pow4(msq)*pow6(mU3)*pow8(mQ3) - 36720*dilog(1 -
         mU32/msq2)*pow4(msq)*pow6(mU3)*pow8(mQ3) +
         11520*logQ3U3*pow4(msq)*pow6(mU3)*pow8(mQ3) - 360*dilog(1 -
         msq2/mU32)* logQ3U3*pow4(msq)*pow6(mU3)*pow8(mQ3) + 19800*dilog(1 -
         mU32/ msq2)*logQ3U3*pow4(msq)*pow6(mU3)*pow8(mQ3) -
         9360*logQ3mq*logQ3U3*pow4(msq)*pow6(mU3)*pow8(mQ3) -
         18360*pow2(logQ3mq)*pow4(msq)*pow6(mU3)*pow8(mQ3) + 9900*logQ3U3*pow2(
         logQ3mq)*pow4(msq)*pow6(mU3)*pow8(mQ3) -
         1080*pow2(logQ3U3)*pow4(msq)*pow6(mU3)*pow8(mQ3) +
         360*logQ3mq*pow2(logQ3U3)*pow4(msq)*pow6(mU3)*pow8(mQ3) -
         180*pow3(logQ3U3)* pow4(msq)*pow6(mU3)*pow8(mQ3) +
         34560*pow4(msq)*pow6(mQ3)*pow8(mU3) - 3600*dilog(1 -
         msq2/mU32)*pow4(msq)*pow6(mQ3)*pow8(mU3) + 48240*dilog(1 -
         mU32/msq2)*pow4(msq)*pow6(mQ3)*pow8(mU3) - 34560*logQ3mq*pow4(
         msq)*pow6(mQ3)*pow8(mU3) - 5760*logQ3U3*pow4(msq)*pow6(mQ3)* pow8(mU3)
         - 360*dilog(1 - msq2/mU32)*logQ3U3*pow4(msq)*pow6(mQ3) *pow8(mU3) +
         14040*dilog(1 - mU32/msq2)*logQ3U3*pow4(msq)*pow6( mQ3)*pow8(mU3) +
         15120*logQ3mq*logQ3U3*pow4(msq)*pow6( mQ3)*pow8(mU3) +
         24120*pow2(logQ3mq)*pow4(msq)*pow6(mQ3)*pow8( mU3) +
         7020*logQ3U3*pow2(logQ3mq)*pow4(msq)*pow6(mQ3)* pow8(mU3) -
         10440*pow2(logQ3U3)*pow4(msq)*pow6(mQ3)*pow8(mU3) +
         360*logQ3mq*pow2(logQ3U3)*pow4(msq)*pow6(mQ3)*pow8(mU3) -
         180*pow3(logQ3U3)*pow4(msq)*pow6(mQ3)*pow8(mU3) + 1800*dilog(1 -
         msq2/mU32)*pow4(mQ3)*pow6(msq)*pow8(mU3) - 27000*dilog(1 - mU32/msq2)*
         pow4(mQ3)*pow6(msq)*pow8(mU3) + 540*dilog(1 - msq2/mU32)*logQ3U3
         *pow4(mQ3)*pow6(msq)*pow8(mU3) + 3420*dilog(1 -
         mU32/msq2)*logQ3U3*pow4(mQ3)*pow6(msq)*pow8(mU3) -
         1800*logQ3mq*logQ3U3*pow4(mQ3)*pow6(msq)*pow8(mU3) -
         13500*pow2(logQ3mq)*pow4( mQ3)*pow6(msq)*pow8(mU3) +
         1710*logQ3U3*pow2(logQ3mq)* pow4(mQ3)*pow6(msq)*pow8(mU3) +
         900*pow2(logQ3U3)*pow4(mQ3)* pow6(msq)*pow8(mU3) -
         540*logQ3mq*pow2(logQ3U3)*pow4(mQ3) *pow6(msq)*pow8(mU3) +
         270*pow3(logQ3U3)*pow4(mQ3)*pow6(msq)* pow8(mU3) -
         1433792*msq2*pow8(mQ3)*pow8(mU3) - 812288*lmQ3MR*msq2*pow8(
         mQ3)*pow8(mU3) - 18000*msq2*dilog(1 - msq2/mU32)*pow8(mQ3)*pow8(mU3) -
         250452*msq2*dilog(1 - mU32/mQ32)*pow8(mQ3)*pow8(mU3) + 107520*lmQ3MR*
         msq2*dilog(1 - mU32/mQ32)*pow8(mQ3)*pow8(mU3) + 102960*msq2*dilog(1 -
         mU32/msq2)*pow8(mQ3)*pow8(mU3) + 353280*msq2*logQ3mq*pow8(mQ3)*
         pow8(mU3) - 9600*lmQ3MR*msq2*logQ3mq*pow8(mQ3)*pow8(mU3) -
         19200*msq2*dilog(1 - mU32/mQ32)*logQ3mq*pow8(mQ3)*pow8(mU3) +
         1226272*msq2*logQ3U3*pow8(mQ3)*pow8(mU3) - 969088*lmQ3MR*msq2*
         logQ3U3*pow8(mQ3)*pow8(mU3) + 1800*msq2*dilog(1 - msq2/mU32)*
         logQ3U3*pow8(mQ3)*pow8(mU3) + 115674*msq2*dilog(1 - mU32/mQ32)*
         logQ3U3*pow8(mQ3)*pow8(mU3) - 35640*msq2*dilog(1 - mU32/msq2)*
         logQ3U3*pow8(mQ3)*pow8(mU3) -
         21360*msq2*logQ3mq*logQ3U3*pow8(mQ3)*pow8(mU3) -
         15360*lmQ3MR*msq2*logQ3mq*logQ3U3*pow8(mQ3)*pow8(mU3) -
         3072*msq2*logU3Q3*pow8(mQ3)* pow8(mU3) -
         384*msq2*logQ3U3*logU3Q3*pow8(mQ3)*pow8(mU3) +
         533760*msq2*pow2(lmQ3MR)*pow8(mQ3)*pow8(mU3) +
         38016*msq2*logQ3U3*pow2(lmQ3MR)*pow8(mQ3)*pow8(mU3) -
         6120*msq2*pow2(logQ3mq) *pow8(mQ3)*pow8(mU3) -
         12060*msq2*logQ3U3*pow2(logQ3mq)* pow8(mQ3)*pow8(mU3) +
         363790*msq2*pow2(logQ3U3)*pow8(mQ3)*pow8( mU3) +
         63552*lmQ3MR*msq2*pow2(logQ3U3)*pow8(mQ3)*pow8(mU3) +
         15480*msq2*logQ3mq*pow2(logQ3U3)*pow8(mQ3)*pow8(mU3) -
         58287*msq2*pow3(logQ3U3)*pow8(mQ3)*pow8(mU3) - 5760*pow4(msq)*
         pow4(mU3)*power10(mQ3) - 1800*dilog(1 -
         msq2/mU32)*pow4(msq)*pow4(mU3)* power10(mQ3) - 1800*dilog(1 -
         mU32/msq2)*pow4(msq)*pow4(mU3)*power10( mQ3) +
         5760*logQ3mq*pow4(msq)*pow4(mU3)*power10(mQ3) +
         2880*logQ3U3*pow4(msq)*pow4(mU3)*power10(mQ3) + 540*dilog(1 -
         msq2/mU32)* logQ3U3*pow4(msq)*pow4(mU3)*power10(mQ3) + 540*dilog(1 -
         mU32/ msq2)*logQ3U3*pow4(msq)*pow4(mU3)*power10(mQ3) -
         1080*logQ3mq*logQ3U3*pow4(msq)*pow4(mU3)*power10(mQ3) -
         900*pow2(logQ3mq)*pow4(msq)*pow4(mU3)*power10(mQ3) + 270*logQ3U3*pow2(
         logQ3mq)*pow4(msq)*pow4(mU3)*power10(mQ3) -
         900*pow2(logQ3U3)*pow4(msq)*pow4(mU3)*power10(mQ3) -
         540*logQ3mq*pow2(logQ3U3)*pow4(msq)*pow4(mU3)*power10(mQ3) +
         270*pow3(logQ3U3) *pow4(msq)*pow4(mU3)*power10(mQ3) -
         1800*mU32*dilog(1 - msq2/mU32)* pow6(msq)*power10(mQ3) -
         1800*mU32*dilog(1 - mU32/msq2)*pow6(msq)* power10(mQ3) +
         540*mU32*dilog(1 - msq2/mU32)*logQ3U3*pow6(msq)* power10(mQ3) +
         540*mU32*dilog(1 - mU32/msq2)*logQ3U3*pow6(msq)* power10(mQ3) +
         1800*mU32*logQ3mq*logQ3U3*pow6(msq)* power10(mQ3) -
         900*mU32*pow2(logQ3mq)*pow6(msq)*power10(mQ3) +
         270*mU32*logQ3U3*pow2(logQ3mq)*pow6(msq)*power10(mQ3) -
         900*mU32*pow2(logQ3U3)*pow6(msq)*power10(mQ3) -
         540*mU32*logQ3mq*pow2(logQ3U3)*pow6(msq)*power10(mQ3) + 270*mU32*pow3(
         logQ3U3)*pow6(msq)*power10(mQ3) + 515872*msq2*pow6(mU3)*power10( mQ3)
         + 499648*lmQ3MR*msq2*pow6(mU3)*power10(mQ3) + 9000*msq2*dilog(1 -
         msq2/mU32)*pow6(mU3)*power10(mQ3) + 42564*msq2*dilog(1 - mU32/mQ32)*
         pow6(mU3)*power10(mQ3) - 53760*lmQ3MR*msq2*dilog(1 - mU32/mQ32)*pow6(
         mU3)*power10(mQ3) - 19800*msq2*dilog(1 - mU32/msq2)*pow6(mU3)*power10(
         mQ3) - 183360*msq2*logQ3mq*pow6(mU3)*power10(mQ3) + 7680*lmQ3MR*
         msq2*logQ3mq*pow6(mU3)*power10(mQ3) + 9600*msq2*dilog(1 - mU32/
         mQ32)*logQ3mq*pow6(mU3)*power10(mQ3) -
         619024*msq2*logQ3U3*pow6(mU3)*power10(mQ3) +
         254688*lmQ3MR*msq2*logQ3U3*pow6( mU3)*power10(mQ3) - 2700*msq2*dilog(1
         - msq2/mU32)*logQ3U3*pow6( mU3)*power10(mQ3) - 28590*msq2*dilog(1 -
         mU32/mQ32)*logQ3U3* pow6(mU3)*power10(mQ3) + 5940*msq2*dilog(1 -
         mU32/msq2)*logQ3U3* pow6(mU3)*power10(mQ3) +
         45240*msq2*logQ3mq*logQ3U3*pow6( mU3)*power10(mQ3) +
         3840*lmQ3MR*msq2*logQ3mq*logQ3U3* pow6(mU3)*power10(mQ3) +
         2688*msq2*logU3Q3*pow6(mU3)*power10( mQ3) -
         576*msq2*logQ3U3*logU3Q3*pow6(mU3)*power10(mQ3) -
         283008*msq2*pow2(lmQ3MR)*pow6(mU3)*power10(mQ3) + 50496*msq2*logQ3U3*
         pow2(lmQ3MR)*pow6(mU3)*power10(mQ3) + 18900*msq2*pow2(logQ3mq)*
         pow6(mU3)*power10(mQ3) - 5670*msq2*logQ3U3*pow2(logQ3mq)*pow6(mU3)*
         power10(mQ3) - 1626*msq2*pow2(logQ3U3)*pow6( mU3)*power10(mQ3) -
         100800*lmQ3MR*msq2*pow2(logQ3U3)*pow6(mU3)* power10(mQ3) -
         1140*msq2*logQ3mq*pow2(logQ3U3)*pow6(mU3)* power10(mQ3) +
         67587*msq2*pow3(logQ3U3)*pow6(mU3)*power10(mQ3) -
         510112*pow8(mU3)*power10(mQ3) - 499648*lmQ3MR*pow8(mU3)*power10(mQ3) -
         5400*dilog(1 - msq2/mU32)*pow8(mU3)*power10(mQ3) - 42564*dilog(1 -
         mU32/mQ32)*pow8(mU3)*power10(mQ3) + 53760*lmQ3MR*dilog(1 - mU32/mQ32)*
         pow8(mU3)*power10(mQ3) + 23400*dilog(1 - mU32/msq2)*pow8(mU3)*power10(
         mQ3) + 177600*logQ3mq*pow8(mU3)*power10(mQ3) - 7680*lmQ3MR*logQ3mq*
         pow8(mU3)*power10(mQ3) - 9600*dilog(1 -
         mU32/mQ32)*logQ3mq*pow8(mU3)*power10(mQ3) +
         616144*logQ3U3*pow8(mU3)*power10( mQ3) -
         254688*lmQ3MR*logQ3U3*pow8(mU3)*power10(mQ3) + 1620* dilog(1 -
         msq2/mU32)*logQ3U3*pow8(mU3)*power10(mQ3) + 28590* dilog(1 -
         mU32/mQ32)*logQ3U3*pow8(mU3)*power10(mQ3) - 7020* dilog(1 -
         mU32/msq2)*logQ3U3*pow8(mU3)*power10(mQ3) - 45960*logQ3mq*
         logQ3U3*pow8(mU3)*power10(mQ3) -
         3840*lmQ3MR*logQ3mq*logQ3U3*pow8(mU3)*power10(mQ3) - 2688*logU3Q3*
         pow8(mU3)*power10(mQ3) + 576*logQ3U3*logU3Q3*pow8(mU3)* power10(mQ3) -
         19200*logU3mq*pow8(mU3)*power10(mQ3) + 5760*logQ3U3*
         logU3mq*pow8(mU3)*power10(mQ3) + 283008*pow2(lmQ3MR)*
         pow8(mU3)*power10(mQ3) - 50496*logQ3U3*pow2(lmQ3MR)*pow8(mU3)*
         power10(mQ3) - 17100*pow2(logQ3mq)*pow8(mU3)*power10(mQ3) +
         5130*logQ3U3*pow2(logQ3mq)*pow8(mU3)*power10(mQ3) + 3426*
         pow2(logQ3U3)*pow8(mU3)*power10(mQ3) + 100800*lmQ3MR*pow2(logQ3U3)*
         pow8(mU3)*power10(mQ3) +
         2220*logQ3mq*pow2(logQ3U3)*pow8(mU3)*power10(mQ3) -
         68127*pow3(logQ3U3)*pow8(mU3)* power10(mQ3) -
         46080*pow4(mQ3)*pow4(msq)*power10(mU3) + 1800*dilog(1 -
         msq2/mU32)*pow4(mQ3)*pow4(msq)*power10(mU3) + 19080*dilog(1 - mU32/
         msq2)*pow4(mQ3)*pow4(msq)*power10(mU3) + 46080*logQ3mq*pow4(mQ3)
         *pow4(msq)*power10(mU3) - 34560*logQ3U3*pow4(mQ3)*pow4(msq)*
         power10(mU3) + 540*dilog(1 - msq2/mU32)*logQ3U3*pow4(mQ3)*pow4(
         msq)*power10(mU3) - 19620*dilog(1 - mU32/msq2)*logQ3U3*pow4(mQ3)
         *pow4(msq)*power10(mU3) + 3960*logQ3mq*logQ3U3*pow4(mQ3)*
         pow4(msq)*power10(mU3) + 9540*pow2(logQ3mq)*pow4(mQ3)*pow4(msq)*
         power10(mU3) - 9810*logQ3U3*pow2(logQ3mq)*pow4(mQ3)*pow4(
         msq)*power10(mU3) + 3780*pow2(logQ3U3)*pow4(mQ3)*pow4(msq)*
         power10(mU3) - 540*logQ3mq*pow2(logQ3U3)*pow4(mQ3)*pow4(
         msq)*power10(mU3) + 270*pow3(logQ3U3)*pow4(mQ3)*pow4(msq)*
         power10(mU3) + 1805120*msq2*pow6(mQ3)*power10(mU3) + 603008*lmQ3MR*
         msq2*pow6(mQ3)*power10(mU3) + 18000*msq2*dilog(1 -
         msq2/mU32)*pow6(mQ3) *power10(mU3) + 425916*msq2*dilog(1 -
         mU32/mQ32)*pow6(mQ3)*power10(mU3) - 107520*lmQ3MR*msq2*dilog(1 -
         mU32/mQ32)*pow6(mQ3)*power10(mU3) - 160560*msq2*dilog(1 -
         mU32/msq2)*pow6(mQ3)*power10(mU3) - 312960*msq2*
         logQ3mq*pow6(mQ3)*power10(mU3) + 19200*msq2*dilog(1 - mU32/mQ32)
         *logQ3mq*pow6(mQ3)*power10(mU3) - 818208*msq2*logQ3U3*
         pow6(mQ3)*power10(mU3) + 1234112*lmQ3MR*msq2*logQ3U3*pow6(mQ3)*
         power10(mU3) + 1800*msq2*dilog(1 - msq2/mU32)*logQ3U3*pow6(mQ3)*
         power10(mU3) - 91098*msq2*dilog(1 - mU32/mQ32)*logQ3U3*pow6(mQ3)
         *power10(mU3) - 6840*msq2*dilog(1 - mU32/msq2)*logQ3U3*pow6(mQ3)
         *power10(mU3) - 36240*msq2*logQ3mq*logQ3U3*pow6(mQ3)* power10(mU3) +
         23040*lmQ3MR*msq2*logQ3mq*logQ3U3*pow6( mQ3)*power10(mU3) +
         768*msq2*logU3Q3*pow6(mQ3)*power10(mU3) +
         768*msq2*logQ3U3*logU3Q3*pow6(mQ3)*power10(mU3) - 480000*
         msq2*pow2(lmQ3MR)*pow6(mQ3)*power10(mU3) - 177024*msq2*logQ3U3*
         pow2(lmQ3MR)*pow6(mQ3)*power10(mU3) - 22680*msq2*pow2(logQ3mq)*
         pow6(mQ3)*power10(mU3) + 2340*msq2*logQ3U3*pow2(logQ3mq)*
         pow6(mQ3)*power10(mU3) - 646682*msq2*pow2(logQ3U3)*pow6(mQ3)*
         power10(mU3) + 206016*lmQ3MR*msq2*pow2(logQ3U3)*pow6(mQ3)*
         power10(mU3) - 11400*msq2*logQ3mq*pow2(logQ3U3)*pow6(mQ3)
         *power10(mU3) - 125385*msq2*pow3(logQ3U3)*pow6(mQ3)*power10(mU3) -
         360*mQ32*dilog(1 - msq2/mU32)*pow6(msq)*power10(mU3) + 16920*mQ32*
         dilog(1 - mU32/msq2)*pow6(msq)*power10(mU3) - 180*mQ32*dilog(1 - msq2/
         mU32)*logQ3U3*pow6(msq)*power10(mU3) + 8460*mQ32*dilog(1 - mU32/
         msq2)*logQ3U3*pow6(msq)*power10(mU3) + 360*mQ32*logQ3mq*
         logQ3U3*pow6(msq)*power10(mU3) + 8460*mQ32*pow2(logQ3mq)*
         pow6(msq)*power10(mU3) + 4230*mQ32*logQ3U3*pow2(logQ3mq)*
         pow6(msq)*power10(mU3) - 180*mQ32*pow2(logQ3U3)*pow6(msq)*
         power10(mU3) + 180*mQ32*logQ3mq*pow2(logQ3U3)*pow6(msq)* power10(mU3)
         - 90*mQ32*pow3(logQ3U3)*pow6(msq)*power10(mU3) +
         1433792*pow8(mQ3)*power10(mU3) + 812288*lmQ3MR*pow8(mQ3)*power10(mU3)
         + 10800*dilog(1 - msq2/mU32)*pow8(mQ3)*power10(mU3) + 250452*dilog(1 -
         mU32/mQ32)*pow8(mQ3)*power10(mU3) - 107520*lmQ3MR*dilog(1 -
         mU32/mQ32)* pow8(mQ3)*power10(mU3) - 75600*dilog(1 -
         mU32/msq2)*pow8(mQ3)*power10( mU3) -
         353280*logQ3mq*pow8(mQ3)*power10(mU3) + 9600*lmQ3MR*logQ3mq*pow8(mQ3)*
         power10(mU3) + 19200*dilog(1 -
         mU32/mQ32)*logQ3mq*pow8(mQ3)*power10(mU3) - 1237792*logQ3U3*pow8(mQ3)*
         power10(mU3) + 969088*lmQ3MR*logQ3U3*pow8(mQ3)*power10(mU3) -
         1080*dilog(1 - msq2/mU32)*logQ3U3*pow8(mQ3)*power10(mU3) -
         115674*dilog(1 - mU32/mQ32)*logQ3U3*pow8(mQ3)*power10(mU3) +
         19080*dilog(1 - mU32/msq2)*logQ3U3*pow8(mQ3)*power10(mU3) +
         34320*logQ3mq*logQ3U3*pow8(mQ3)*power10(mU3) + 15360*
         lmQ3MR*logQ3mq*logQ3U3*pow8(mQ3)*power10(mU3) +
         3072*logU3Q3*pow8(mQ3)*power10(mU3) + 384*logQ3U3*logU3Q3*
         pow8(mQ3)*power10(mU3) + 38400*logU3mq*pow8(mQ3)*power10(mU3) -
         3840*logQ3U3*logU3mq*pow8(mQ3)*power10(mU3) - 533760*
         pow2(lmQ3MR)*pow8(mQ3)*power10(mU3) - 38016*logQ3U3*pow2(lmQ3MR)
         *pow8(mQ3)*power10(mU3) + 19800*pow2(logQ3mq)*pow8(mQ3)*power10( mU3)
         + 3780*logQ3U3*pow2(logQ3mq)*pow8(mQ3)*power10(mU3) -
         364510*pow2(logQ3U3)*pow8(mQ3)*power10(mU3) - 63552*lmQ3MR*
         pow2(logQ3U3)*pow8(mQ3)*power10(mU3) - 16200*logQ3mq*
         pow2(logQ3U3)*pow8(mQ3)*power10(mU3) + 58647*pow3(logQ3U3
         )*pow8(mQ3)*power10(mU3)))/(36.*mQ32*(-msq2 +
         mU32)*pow4(mU3)*pow7(mQ32 - mU32))+ (2*Xt5*(88*mU32*dilog(1 -
         mU32/mQ32)*pow12(mQ3) + 48*mU32*logQ3U3* pow12(mQ3) - 32*mU32*dilog(1
         - mU32/mQ32)*logQ3U3*pow12(mQ3) + 48*mU32*logU3Q3*pow12(mQ3) -
         24*mU32*logQ3U3*logU3Q3*pow12(mQ3) + 6848*mQ32*pow12(mU3) -
         4416*lmQ3MR*mQ32*pow12(mU3) - 120*mQ32*dilog(1 - msq2/mU32)*pow12(mU3)
         + 612*mQ32*dilog(1 - mU32/ mQ32)*pow12(mU3) - 120*mQ32*dilog(1 -
         mU32/msq2)*pow12(mU3) - 160*mQ32* logQ3mq*pow12(mU3) +
         11408*mQ32*logQ3U3*pow12(mU3) - 5536*lmQ3MR*mQ32*logQ3U3*pow12(mU3) -
         60*mQ32*dilog(1 - msq2/ mU32)*logQ3U3*pow12(mU3) + 562*mQ32*dilog(1 -
         mU32/mQ32)*logQ3U3*pow12(mU3) - 60*mQ32*dilog(1 - mU32/msq2)*logQ3U3*
         pow12(mU3) + 40*mQ32*logQ3mq*logQ3U3*pow12(mU3) + 112*
         mQ32*logU3Q3*pow12(mU3) + 56*mQ32*logQ3U3*logU3Q3* pow12(mU3) +
         480*mQ32*logU3mq*pow12(mU3) + 240*mQ32*logQ3U3*logU3mq*pow12(mU3) -
         12*dilog(1 - mU32/mQ32)*pow14(mQ3) + 6*dilog(1 -
         mU32/mQ32)*logQ3U3*pow14(mQ3) + 512*pow14(mU3) - 512*lmQ3MR*pow14(mU3)
         + 512*logQ3U3*pow14(mU3) - 60*mQ32*pow12( mU3)*pow2(logQ3mq) -
         30*mQ32*logQ3U3*pow12(mU3)*pow2(logQ3mq) +
         20*mU32*pow12(mQ3)*pow2(logQ3U3) + 6158*mQ32* pow12(mU3)*pow2(logQ3U3)
         - 512*lmQ3MR*mQ32*pow12(mU3)*pow2(logQ3U3) +
         60*mQ32*logQ3mq*pow12(mU3)*pow2(logQ3U3) - 6*pow14(mQ3)*pow2(logQ3U3)
         - 16*mU32*pow12(mQ3)*pow3(logQ3U3) + 763*mQ32*pow12(mU3)*pow3(logQ3U3)
         + 3*pow14(mQ3)*pow3( logQ3U3) + 2*mQ32*dilog(1 - mQ32/mU32)*(-2*mQ32 +
         2*mU32 + (mQ32 + mU32)*logQ3U3)*pow3(mQ32 - mU32)*(-10*mQ32*mU32 +
         3*pow4(mQ3) - 57*pow4(mU3)) + 960*(mQ32 - mU32)*dilog(1 -
         msq2/mQ32)*(-2*mQ32 + 2* mU32 + (mQ32 +
         mU32)*logQ3U3)*pow2(-(mQ3*msq2) + pow3(mQ3))* pow4(mU3) - 720*dilog(1
         - msq2/mU32)*pow4(msq)*pow4(mU3)*pow6(mQ3) - 2640*dilog(1 -
         mU32/msq2)*pow4(msq)*pow4(mU3)*pow6(mQ3) + 960*dilog(1 -
         mU32/msq2)*logQ3U3*pow4(msq)*pow4(mU3)*pow6(mQ3) +
         720*logQ3mq*logQ3U3*pow4(msq)*pow4(mU3)*pow6(mQ3) -
         1320*pow2(logQ3mq)*pow4(msq)*pow4(mU3)*pow6(mQ3) + 480*logQ3U3*pow2(
         logQ3mq)*pow4(msq)*pow4(mU3)*pow6(mQ3) - 360*pow2(logQ3U3
         )*pow4(msq)*pow4(mU3)*pow6(mQ3) + 480*dilog(1 - msq2/mU32)*pow4(mQ3)*
         pow4(msq)*pow6(mU3) + 4320*dilog(1 - mU32/msq2)*pow4(mQ3)*pow4(msq)*
         pow6(mU3) + 120*dilog(1 - msq2/mU32)*logQ3U3*pow4(mQ3)*pow4(msq)
         *pow6(mU3) + 120*dilog(1 - mU32/msq2)*logQ3U3*pow4(mQ3)*pow4(
         msq)*pow6(mU3) - 480*logQ3mq*logQ3U3*pow4(mQ3)*pow4(msq)* pow6(mU3) +
         2160*pow2(logQ3mq)*pow4(mQ3)*pow4(msq)*pow6(mU3) +
         60*logQ3U3*pow2(logQ3mq)*pow4(mQ3)*pow4(msq)*pow6(mU3) +
         240*pow2(logQ3U3)*pow4(mQ3)*pow4(msq)*pow6(mU3) -
         120*logQ3mq*pow2(logQ3U3)*pow4(mQ3)*pow4(msq)*pow6(mU3) +
         60*pow3(logQ3U3)*pow4(mQ3)*pow4(msq)*pow6(mU3) -
         5760*msq2*pow6(mQ3)*pow6( mU3) + 1440*msq2*dilog(1 -
         msq2/mU32)*pow6(mQ3)*pow6(mU3) - 6240*msq2* dilog(1 -
         mU32/msq2)*pow6(mQ3)*pow6(mU3) + 2880*msq2*logQ3mq*
         pow6(mQ3)*pow6(mU3) + 1920*msq2*logQ3U3*pow6(mQ3)*pow6(mU3) -
         2400*msq2*logQ3mq*logQ3U3*pow6(mQ3)*pow6(mU3) + 2880*
         msq2*logU3mq*pow6(mQ3)*pow6(mU3) - 3120*msq2*pow2(logQ3mq
         )*pow6(mQ3)*pow6(mU3) + 1680*msq2*pow2(logQ3U3)*pow6(mQ3)*pow6( mU3) +
         480*mU32*dilog(1 - msq2/mU32)*pow4(msq)*pow8(mQ3) + 480*mU32* dilog(1
         - mU32/msq2)*pow4(msq)*pow8(mQ3) - 120*mU32*dilog(1 - msq2/
         mU32)*logQ3U3*pow4(msq)*pow8(mQ3) - 120*mU32*dilog(1 - mU32/
         msq2)*logQ3U3*pow4(msq)*pow8(mQ3) -
         480*mU32*logQ3mq*logQ3U3*pow4(msq)*pow8(mQ3) +
         240*mU32*pow2(logQ3mq)*pow4( msq)*pow8(mQ3) -
         60*mU32*logQ3U3*pow2(logQ3mq)*pow4(msq)* pow8(mQ3) +
         240*mU32*pow2(logQ3U3)*pow4(msq)*pow8(mQ3) + 120*
         mU32*logQ3mq*pow2(logQ3U3)*pow4(msq)*pow8(mQ3) - 60*mU32*
         pow3(logQ3U3)*pow4(msq)*pow8(mQ3) + 1920*msq2*pow4(mU3)*pow8( mQ3) -
         960*msq2*dilog(1 - msq2/mU32)*pow4(mU3)*pow8(mQ3) + 2880*msq2* dilog(1
         - mU32/msq2)*pow4(mU3)*pow8(mQ3) - 2880*msq2*logQ3U3*
         pow4(mU3)*pow8(mQ3) + 240*msq2*dilog(1 - msq2/mU32)*logQ3U3*
         pow4(mU3)*pow8(mQ3) - 1680*msq2*dilog(1 - mU32/msq2)*logQ3U3*
         pow4(mU3)*pow8(mQ3) + 1440*msq2*logQ3mq*logQ3U3*pow4(mU3) *pow8(mQ3) -
         1920*msq2*logU3mq*pow4(mU3)*pow8(mQ3) + 480*msq2*
         logQ3U3*logU3mq*pow4(mU3)*pow8(mQ3) +
         1440*msq2*pow2(logQ3mq)*pow4(mU3)*pow8(mQ3) -
         840*msq2*logQ3U3*pow2(logQ3mq)*pow4(mU3)*pow8(mQ3) -
         240*msq2*logQ3mq*pow2(logQ3U3)*pow4(mU3)*pow8(mQ3) +
         120*msq2*pow3(logQ3U3)*pow4( mU3)*pow8(mQ3) -
         32000*pow6(mU3)*pow8(mQ3) + 25728*lmQ3MR*pow6(mU3)* pow8(mQ3) +
         480*dilog(1 - msq2/mU32)*pow6(mU3)*pow8(mQ3) - 2160*dilog(1 -
         mU32/mQ32)*pow6(mU3)*pow8(mQ3) + 4320*dilog(1 - mU32/msq2)*pow6(mU3)*
         pow8(mQ3) + 640*logQ3mq*pow6(mU3)*pow8(mQ3) +
         8672*logQ3U3*pow6(mU3)*pow8(mQ3) -
         19840*lmQ3MR*logQ3U3*pow6(mU3)*pow8( mQ3) - 120*dilog(1 -
         msq2/mU32)*logQ3U3*pow6(mU3)*pow8(mQ3) - 400*dilog(1 -
         mU32/mQ32)*logQ3U3*pow6(mU3)*pow8(mQ3) - 120* dilog(1 -
         mU32/msq2)*logQ3U3*pow6(mU3)*pow8(mQ3) -
         3520*logQ3mq*logQ3U3*pow6(mU3)*pow8(mQ3) - 160*logU3Q3*
         pow6(mU3)*pow8(mQ3) + 112*logQ3U3*logU3Q3*pow6(mU3)*pow8( mQ3) -
         1920*logU3mq*pow6(mU3)*pow8(mQ3) + 480*logQ3U3*
         logU3mq*pow6(mU3)*pow8(mQ3) + 2160*pow2(logQ3mq)*pow6( mU3)*pow8(mQ3)
         - 60*logQ3U3*pow2(logQ3mq)*pow6(mU3)*pow8( mQ3) +
         9256*pow2(logQ3U3)*pow6(mU3)*pow8(mQ3) + 896*lmQ3MR*pow2(
         logQ3U3)*pow6(mU3)*pow8(mQ3) -
         40*logQ3mq*pow2(logQ3U3)*pow6(mU3)*pow8(mQ3) -
         212*pow3(logQ3U3)*pow6(mU3)*pow8( mQ3) + 5760*msq2*pow4(mQ3)*pow8(mU3)
         - 960*msq2*dilog(1 - msq2/mU32)* pow4(mQ3)*pow8(mU3) +
         2880*msq2*dilog(1 - mU32/msq2)*pow4(mQ3)*pow8( mU3) -
         3840*msq2*logQ3mq*pow4(mQ3)*pow8(mU3) +
         2880*msq2*logQ3U3*pow4(mQ3)*pow8(mU3) - 240*msq2*dilog(1 -
         msq2/mU32)*logQ3U3*pow4(mQ3)*pow8(mU3) + 1680*msq2*dilog(1 -
         mU32/msq2)*logQ3U3*pow4(mQ3)*pow8(mU3) + 480*msq2*logQ3mq*logQ3U3
         *pow4(mQ3)*pow8(mU3) - 1920*msq2*logU3mq*pow4(mQ3)*pow8(mU3) -
         480*msq2*logQ3U3*logU3mq*pow4(mQ3)*pow8(mU3) + 1440*msq2*
         pow2(logQ3mq)*pow4(mQ3)*pow8(mU3) + 840*msq2*logQ3U3*
         pow2(logQ3mq)*pow4(mQ3)*pow8(mU3) - 960*msq2*pow2(logQ3U3
         )*pow4(mQ3)*pow8(mU3) + 240*msq2*logQ3mq*pow2(logQ3U3)*
         pow4(mQ3)*pow8(mU3) - 120*msq2*pow3(logQ3U3)*pow4(mQ3)*pow8(mU3) -
         120*mQ32*dilog(1 - msq2/mU32)*pow4(msq)*pow8(mU3) - 2040*mQ32*dilog(1
         - mU32/msq2)*pow4(msq)*pow8(mU3) - 60*mQ32*dilog(1 -
         msq2/mU32)*logQ3U3*pow4(msq)*pow8(mU3) - 1020*mQ32*dilog(1 -
         mU32/msq2)*logQ3U3*pow4(msq)*pow8(mU3) + 120*mQ32*logQ3mq*logQ3U3
         *pow4(msq)*pow8(mU3) - 1020*mQ32*pow2(logQ3mq)*pow4(msq)*pow8( mU3) -
         510*mQ32*logQ3U3*pow2(logQ3mq)*pow4(msq)*pow8(mU3) -
         60*mQ32*pow2(logQ3U3)*pow4(msq)*pow8(mU3) +
         60*mQ32*logQ3mq*pow2(logQ3U3)*pow4(msq)*pow8(mU3) -
         30*mQ32*pow3(logQ3U3)*pow4(msq)*pow8(mU3) + 48256*pow6(mQ3)*pow8(mU3)
         - 37120*lmQ3MR* pow6(mQ3)*pow8(mU3) - 720*dilog(1 -
         msq2/mU32)*pow6(mQ3)*pow8(mU3) + 3500*dilog(1 -
         mU32/mQ32)*pow6(mQ3)*pow8(mU3) - 2640*dilog(1 - mU32/
         msq2)*pow6(mQ3)*pow8(mU3) - 960*logQ3mq*pow6(mQ3)*pow8(mU3) +
         14240*logQ3U3*pow6(mQ3)*pow8(mU3) + 8512*lmQ3MR*logQ3U3*
         pow6(mQ3)*pow8(mU3) + 1490*dilog(1 - mU32/mQ32)*logQ3U3*pow6(
         mQ3)*pow8(mU3) - 960*dilog(1 - mU32/msq2)*logQ3U3*pow6(mQ3)* pow8(mU3)
         + 1680*logQ3mq*logQ3U3*pow6(mQ3)*pow8(mU3) +
         480*logU3Q3*pow6(mQ3)*pow8(mU3) -
         48*logQ3U3*logU3Q3*pow6(mQ3)*pow8(mU3) +
         2880*logU3mq*pow6(mQ3)*pow8(mU3) -
         1320*pow2(logQ3mq)*pow6(mQ3)*pow8(mU3) - 480*logQ3U3*
         pow2(logQ3mq)*pow6(mQ3)*pow8(mU3) - 9730*pow2(logQ3U3)*
         pow6(mQ3)*pow8(mU3) + 5056*lmQ3MR*pow2(logQ3U3)*pow6(mQ3)*pow8( mU3) +
         800*logQ3mq*pow2(logQ3U3)*pow6(mQ3)*pow8(mU3) -
         2695*pow3(logQ3U3)*pow6(mQ3)*pow8(mU3) + 240*msq2*mU32*dilog(1 -
         msq2/mU32)*power10(mQ3) + 240*msq2*mU32*dilog(1 - mU32/msq2)*power10(
         mQ3) - 480*msq2*mU32*logQ3mq*power10(mQ3) +
         480*msq2*mU32*logQ3U3*power10(mQ3) - 120*msq2*mU32*dilog(1 -
         msq2/mU32)*logQ3U3*power10(mQ3) - 120*msq2*mU32*dilog(1 -
         mU32/msq2)*logQ3U3* power10(mQ3) + 480*msq2*mU32*logU3mq*power10(mQ3)
         - 240*msq2* mU32*logQ3U3*logU3mq*power10(mQ3) + 120*msq2*mU32*pow2(
         logQ3mq)*power10(mQ3) -
         60*msq2*mU32*logQ3U3*pow2(logQ3mq)*power10(mQ3) -
         120*msq2*mU32*pow2(logQ3U3)*power10( mQ3) +
         120*msq2*mU32*logQ3mq*pow2(logQ3U3)*power10(mQ3) -
         60*msq2*mU32*pow3(logQ3U3)*power10(mQ3) - 120*dilog(1 - msq2/
         mU32)*pow4(msq)*power10(mQ3) - 120*dilog(1 - mU32/msq2)*pow4(msq)*
         power10(mQ3) + 60*dilog(1 - msq2/mU32)*logQ3U3*pow4(msq)* power10(mQ3)
         + 60*dilog(1 - mU32/msq2)*logQ3U3*pow4(msq)* power10(mQ3) +
         120*logQ3mq*logQ3U3*pow4(msq)*power10(mQ3) -
         60*pow2(logQ3mq)*pow4(msq)*power10(mQ3) + 30*logQ3U3*
         pow2(logQ3mq)*pow4(msq)*power10(mQ3) - 60*pow2(logQ3U3)*
         pow4(msq)*power10(mQ3) - 60*logQ3mq*pow2(logQ3U3)*pow4(
         msq)*power10(mQ3) + 30*pow3(logQ3U3)*pow4(msq)*power10(mQ3) +
         7872*pow4(mU3)*power10(mQ3) - 6592*lmQ3MR*pow4(mU3)*power10(mQ3) -
         120* dilog(1 - msq2/mU32)*pow4(mU3)*power10(mQ3) + 380*dilog(1 -
         mU32/mQ32)* pow4(mU3)*power10(mQ3) - 2040*dilog(1 -
         mU32/msq2)*pow4(mU3)*power10( mQ3) -
         160*logQ3mq*pow4(mU3)*power10(mQ3) - 6448*logQ3U3*
         pow4(mU3)*power10(mQ3) + 8544*lmQ3MR*logQ3U3*pow4(mU3)*power10( mQ3) +
         60*dilog(1 - msq2/mU32)*logQ3U3*pow4(mU3)*power10(mQ3) - 10*dilog(1 -
         mU32/mQ32)*logQ3U3*pow4(mU3)*power10(mQ3) + 1020* dilog(1 -
         mU32/msq2)*logQ3U3*pow4(mU3)*power10(mQ3) +
         1800*logQ3mq*logQ3U3*pow4(mU3)*power10(mQ3) - 80*logU3Q3*
         pow4(mU3)*power10(mQ3) - 8*logQ3U3*logU3Q3*pow4(mU3)* power10(mQ3) +
         480*logU3mq*pow4(mU3)*power10(mQ3) -
         240*logQ3U3*logU3mq*pow4(mU3)*power10(mQ3) -
         1020*pow2(logQ3mq)*pow4(mU3)*power10(mQ3) + 510*logQ3U3*pow2(logQ3mq)
         *pow4(mU3)*power10(mQ3) - 1574*pow2(logQ3U3)*pow4(mU3)*power10( mQ3) -
         2496*lmQ3MR*pow2(logQ3U3)*pow4(mU3)*power10(mQ3) - 860*
         logQ3mq*pow2(logQ3U3)*pow4(mU3)*power10(mQ3) + 1161*pow3(
         logQ3U3)*pow4(mU3)*power10(mQ3) - 1920*mQ32*msq2*power10(mU3) +
         240*mQ32*msq2*dilog(1 - msq2/mU32)*power10(mU3) +
         240*mQ32*msq2*dilog(1 - mU32/msq2)*power10(mU3) +
         1440*mQ32*msq2*logQ3mq*power10(mU3) -
         2400*mQ32*msq2*logQ3U3*power10(mU3) + 120*mQ32*msq2*dilog(1 -
         msq2/mU32)*logQ3U3*power10(mU3) + 120*mQ32*msq2*dilog(1 - mU32/
         msq2)*logQ3U3*power10(mU3) +
         480*mQ32*msq2*logQ3mq*logQ3U3*power10(mU3) +
         480*mQ32*msq2*logU3mq*power10(mU3) +
         240*mQ32*msq2*logQ3U3*logU3mq*power10(mU3) + 120*mQ32*
         msq2*pow2(logQ3mq)*power10(mU3) + 60*mQ32*msq2*logQ3U3*
         pow2(logQ3mq)*power10(mU3) - 600*mQ32*msq2*pow2(logQ3U3)* power10(mU3)
         - 120*mQ32*msq2*logQ3mq*pow2(logQ3U3)* power10(mU3) +
         60*mQ32*msq2*pow3(logQ3U3)*power10(mU3) - 31488*
         pow4(mQ3)*power10(mU3) + 22912*lmQ3MR*pow4(mQ3)*power10(mU3) + 480*
         dilog(1 - msq2/mU32)*pow4(mQ3)*power10(mU3) - 2408*dilog(1 -
         mU32/mQ32) *pow4(mQ3)*power10(mU3) + 480*dilog(1 -
         mU32/msq2)*pow4(mQ3)*power10( mU3) +
         640*logQ3mq*pow4(mQ3)*power10(mU3) - 28432*logQ3U3
         *pow4(mQ3)*power10(mU3) + 8320*lmQ3MR*logQ3U3*pow4(mQ3)*power10( mU3)
         + 120*dilog(1 - msq2/mU32)*logQ3U3*pow4(mQ3)*power10(mU3) -
         1616*dilog(1 - mU32/mQ32)*logQ3U3*pow4(mQ3)*power10(mU3) + 120*
         dilog(1 - mU32/msq2)*logQ3U3*pow4(mQ3)*power10(mU3) -
         400*logU3Q3*pow4(mQ3)*power10(mU3) - 88*logQ3U3*logU3Q3*
         pow4(mQ3)*power10(mU3) - 1920*logU3mq*pow4(mQ3)*power10(mU3) -
         480*logQ3U3*logU3mq*pow4(mQ3)*power10(mU3) + 240*pow2(
         logQ3mq)*pow4(mQ3)*power10(mU3) +
         60*logQ3U3*pow2(logQ3mq)*pow4(mQ3)*power10(mU3) -
         4124*pow2(logQ3U3)*pow4( mQ3)*power10(mU3) -
         2944*lmQ3MR*pow2(logQ3U3)*pow4(mQ3)*power10( mU3) +
         40*logQ3mq*pow2(logQ3U3)*pow4(mQ3)*power10(mU3) +
         1508*pow3(logQ3U3)*pow4(mQ3)*power10(mU3)))/(3.*mQ32*pow3(mU3)*
         pow7(mQ32 - mU32))+ (1280*mU32*Xt6*(-2*mQ32 + 2*mU32 + (mQ32 +
         mU32)*logQ3U3)*pow2(- mQ32 + mU32 + mQ32*logQ3U3))/(3.*pow7(mQ32 -
         mU32));
      }

      case Limits::MQ3_EQ_MU3_EQ_M3: {
         const double logQ3mq = std::log(mQ32/msq2);
         const double logmqQ3 = -logQ3mq;

         return
         (Xt4*(10990*msq2*pow4(mQ3) - 296*lmQ3MR*msq2*pow4(mQ3) - 2724*
         msq2*logQ3mq*pow4(mQ3) - 480*lmQ3MR*msq2*logQ3mq*pow4(mQ3) - 1044*
         msq2*logmqQ3*pow4(mQ3) - 312*msq2*pow2(lmQ3MR)*pow4(mQ3) + 612*
         msq2*pow2(logQ3mq)*pow4(mQ3) - 432*mQ32*logQ3mq*pow4(msq) -
         432*mQ32*logmqQ3*pow4(msq) + 279*mQ32*pow2(logQ3mq)* pow4(msq) -
         10990*pow6(mQ3) + 296*lmQ3MR*pow6(mQ3) + 1200*logQ3mq*pow6(mQ3) +
         480*lmQ3MR*logQ3mq*pow6(mQ3) + 312*pow2( lmQ3MR)*pow6(mQ3) -
         81*pow2(logQ3mq)*pow6(mQ3) + 18*dilog(1 -
         mQ32/msq2)*(28*msq2*pow4(mQ3) + 31*mQ32*pow4(msq) + 31*pow6(mQ3) - 90*
         pow6(msq)) + 2676*logQ3mq*pow6(msq) + 2676*logmqQ3*pow6( msq) -
         810*pow2(logQ3mq)*pow6(msq) - 18*dilog(1 - msq2/mQ32)*(-
         68*msq2*pow4(mQ3) - 31*mQ32*pow4(msq) + 9*pow6(mQ3) +
         90*pow6(msq))))/( 27.*(mQ32 - msq2)*pow8(mQ3))
         -(Xt5*(2176*msq2*pow4(mQ3) - 11264*lmQ3MR*msq2*pow4(mQ3) - 3744*msq2*
         logQ3mq*pow4(mQ3) + 2016*msq2*logmqQ3*pow4(mQ3) + 2331*msq2*
         pow2(logQ3mq)*pow4(mQ3) - 47232*mQ32*logQ3mq*pow4(msq) -
         47232*mQ32*logmqQ3*pow4(msq) + 14085*mQ32*pow2(logQ3mq)* pow4(msq) -
         2176*pow6(mQ3) + 11264*lmQ3MR*pow6(mQ3) - 1920*logQ3mq* pow6(mQ3) -
         1449*pow2(logQ3mq)*pow6(mQ3) + 52896*logQ3mq*pow6(msq) +
         52896*logmqQ3*pow6(msq) - 14967*pow2(logQ3mq)* pow6(msq) - 18*dilog(1
         - mQ32/msq2)*(-259*msq2*pow4(mQ3) - 1565* mQ32*pow4(msq) +
         161*pow6(mQ3) + 1663*pow6(msq)) - 18*dilog(1 - msq2/
         mQ32)*(-259*msq2*pow4(mQ3) - 1565*mQ32*pow4(msq) + 161*pow6(mQ3) +
         1663*pow6(msq))))/(216.*(mQ32 - msq2)*pow9(mQ3))+
         (160*Xt6)/(9.*pow6(mQ3));
      }

      case Limits::DEGENERATE: {
         return (2*Xt4*(-5735 + 148*lmQ3MR + 156*pow2(lmQ3MR)))/(27.*pow4(mQ3))+
         (-176*(-7 + 8*lmQ3MR)*Xt5)/(27.*pow5(mQ3))
         +160/(9.*pow6(mQ3))*Xt6;
      }

      default:
         break;
   };

   throw std::runtime_error("Mass limit not included!");
}


/**
 * Sets the mass limit to check terms
 * @param limit an integer key for a mass limit
 */
void ThresholdCalculator::setLimit(Limits limit)
{
   p.massLimit3LThreshold = static_cast<int>(limit);
}

/**
 * Get the mass limit determined by ThresholdCalculator
 * @return The determined mass limit
 */
Limits ThresholdCalculator::getLimit() const
{
   return static_cast<Limits>(p.massLimit3LThreshold);
}

void ThresholdCalculator::setXtOrderOfDeltaLambdaAtAs2(int xtOrder)
{
   xtOrderLambdaAtAs2 = xtOrder;
}

} // namespace mh2_eft
} // namespace himalaya


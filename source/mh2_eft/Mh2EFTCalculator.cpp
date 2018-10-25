// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "Mh2EFTCalculator.hpp"
#include "ThresholdCalculator.hpp"
#include "EFTFlags.hpp"
#include "Logger.hpp"
#include "dilog.h"
#include <cmath>
#include <string>

namespace himalaya {
namespace mh2_eft {

namespace {
   const double zt3 = 1.2020569031595942853997381615114;
   const double Pi  = 3.1415926535897932384626433832795;
   const double log2 = std::log(2.);
   template <typename T> T pow2(T x)  { return x*x; }
   template <typename T> T pow3(T x)  { return x*x*x; }
   template <typename T> T pow4(T x)  { return x*x*x*x; }
   template <typename T> T pow5(T x)  { return x*x*x*x*x; }
   template <typename T> T pow6(T x)  { return x*x*x*x*x*x; }
   template <typename T> T pow7(T x)  { return x*x*x*x*x*x*x; }
   template <typename T> T pow8(T x)  { return x*x*x*x*x*x*x*x; }
   template <typename T> T pow9(T x)  { return x*x*x*x*x*x*x*x*x; }
   template <typename T> T power10(T x) { return pow9(x)*x; }
   template <typename T> T pow11(T x) { return pow2(x)*pow9(x); }
   template <typename T> T pow12(T x) { return pow2(pow6(x)); }
   template <typename T> T pow13(T x) { return pow4(x)*pow9(x); }
   template <typename T> T pow14(T x) { return pow2(pow7(x)); }
   template <typename T> T pow15(T x) { return pow6(x)*pow9(x); }
   template <typename T> T pow16(T x) { return pow2(pow8(x)); }
   template <typename T> T pow17(T x) { return pow16(x)*x; }
   template <typename T> T pow18(T x) { return pow2(pow9(x)); }
   template <typename T> T pow19(T x) { return pow18(x)*x; }
   template <typename T> T pow20(T x) { return pow19(x)*x; }

double isNaN(double var, const std::string& msg = "")
{
   if (std::isnan(var)) {
      WARNING_MSG("NaN appeared in calculation of threshold correction"
                  << (msg.empty() ? "" : " " + msg) << "!");
      return 0.;
   }
   return var;
}

double F5(double x) noexcept {
   if(x == 1.) return 1.;

   const double x2 = pow2(x);

   return 3*x*(1 - pow2(x2) + 2*x2*log(x2))/pow3(1 - x2);
}

double F6(double x) noexcept {
   if(x == 1.) return 0.;

   const double x2 = pow2(x);

   return (x2 - 3)/(4.*(1-x2)) + x2*(x2 - 2)*log(x2)/(2*pow2(1 - x2));
}

double F9(double x1, double x2) noexcept {
   if(x1 == 1. && x2 == 1.){
      return 1.;
   }
   if(x1 == 1. || x2 == 1.){
      if (x1 == 1.) x1 = x2;

      const double x12 = pow2(x1);
      return 2*(1 - x12 + x12*log(x12))/pow2(x12 - 1);
   }
   if(x1 == x2){
      const double x12 = pow2(x1);
      return 2*(-1 + x12 - log(x12))/pow2(x12 - 1);
   }

   const double x12 = pow2(x1);
   const double x22 = pow2(x2);

   return 2*(x12*log(x12)/(x12 - 1) - x22*log(x22)/(x22 - 1))/(x12 - x22);
}

} // anonymous namespace

/**
 *        Constructor
 *         @param p_ a HimalayaInterface struct
 *         @param msq2_ the averaged squark mass of the first two generations squared
 *         @param verbose a bool enable the output of the parameter validation. Enabled by default
 */
Mh2EFTCalculator::Mh2EFTCalculator(
   const himalaya::Parameters& p_, double msq2_, bool verbose)
   : p(p_), msq2(msq2_)
{
   p.validate(verbose);

   if (!std::isfinite(msq2_))
      msq2 = p.calculateMsq2();

   // fill order map
   orderMap.clear();
   for (int i = EFTOrders::FIRST; i < EFTOrders::NUMBER_OF_EFT_ORDERS; i++) {
      orderMap.emplace(i, 1u);
   }
}

void Mh2EFTCalculator::setCorrectionFlag(int variable, int enable)
{
   if(enable < 0 || enable > 1) INFO_MSG("You can only enable (1) or disable (0) corrections!");
   if(orderMap.find(variable) == orderMap.end()) INFO_MSG("Your variable is not defined in the EFTOrders enum!");
   orderMap.at(variable) = enable;
}

double Mh2EFTCalculator::getDeltaMh2EFT0Loop() const
{
   return pow2(p.MZ * std::cos(2 * std::atan(p.vu/p.vd)));
}

/**
 *         Returns the 1-loop EFT contribution to the Higgs mass
 *         @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 *         @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 */
double Mh2EFTCalculator::getDeltaMh2EFT1Loop(int omitSMLogs, int omitMSSMLogs) const
{
   ThresholdCalculator thresholdCalculator(p, msq2);

   using std::log;
   const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));

   const double v2 = pow2(p.vu) + pow2(p.vd);
   const double gt = sqrt(2)*p.Mt/std::sqrt(v2);

   // 1-Loop prefactor at
   const double pref_at = 1./pow2(4*Pi) * pow2(p.Mt * gt);

   const double beta = atan(p.vu/p.vd);
   const double cbeta = cos(beta);
   const double c2beta = cos(2*beta);
   const double sbeta = sin(beta);
   const double mhtree = std::abs(c2beta*p.MZ);
   const double yt = sqrt(2.)*p.Mt/p.vu;
   const double yb = sqrt(2.)*p.Mb/p.vd;
   const double ytau = sqrt(2.)*p.Mtau/p.vd;
   const double lmhtreeMt = log(pow2(mhtree / p.Mt));
   const double lmwMt = log(pow2(p.MW / p.Mt));
   const double lmzMt = log(pow2(p.MZ / p.Mt));
   const int Xi = 1;        // gauge parameter

   // Threshold corrections
   const double dlambdayb2g12 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YB2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdag14 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_G14, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaregg14 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_REG_G14, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdachig14 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_CHI_G14, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dg1g1 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::G1_G1, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdachig24 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_CHI_G24, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdag24 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_G24, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dg2g2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::G2_G2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaregg24 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_REG_G24, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdag12g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_G12_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaregg12g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_REG_G12_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdachig12g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_CHI_G12_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb2g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YB2_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb4 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YB4, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayt2g12 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YT2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dvyt2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::VEV_YT2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayt2g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YT2_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaytau2g12 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaytau2g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU2_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaytau4 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU4, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dvg12 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::VEV_G12, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dvg22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::VEV_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dvyb2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::VEV_YB2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dvytau2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::VEV_YTAU2, RenSchemes::DRBARPRIME, omitMSSMLogs);

   const double bbhDR = 2 - Pi/sqrt(3.) - log(pow2(mhtree/p.Mt));

   const double bbwDR =
      2 - lmwMt - 2*std::asin(mhtree/(2.*p.MW))*sqrt(-1 + 4*pow2(p.MW/mhtree));

   const double bbzDR =
      2 - lmzMt - 2*std::asin(sqrt(pow2(c2beta))/2.)*sqrt(-1 + 4/pow2(c2beta));

   // corrections to Mh2
   const double dmh2g12g22 = isNaN(orderMap.at(EFTOrders::G12G22)*((v2*(120 + 200
      *dlambdag12g22 + 200*dlambdaregg12g22 + 200*dlambdachig12g22 - 180*lmMt
      + 200*dvg12*pow2(c2beta) + 120*dvg22*pow2(c2beta) + 180*lmMt*pow2(c2beta)
      - 60*lmMt*Xi*pow2(c2beta) + 30*(-2 + pow2(c2beta))*pow2(c2beta)*(-2
      + lmwMt + 2*asin(mhtree/(2.*p.MW))*sqrt(-1 + (4*pow2(p.MW))/pow2(mhtree)))
      - 18*pow4(c2beta) - 180*lmMt*pow4(c2beta) + 45*(-6 + 3*lmhtreeMt + sqrt(3.)
      *Pi)*pow4(c2beta) + 15*(-2 + lmzMt + 2*asin(std::abs(c2beta)/2.)*sqrt(-1
      + 4/pow2(c2beta)))*(12 - 4* pow2(c2beta) + pow4(c2beta))))/400.),
      "dmh2g12g22");
   const double dmh2g14 = isNaN(orderMap.at(EFTOrders::G14)*(-(v2*(-20*(18 + 100
      *dlambdag14 + 100*dlambdaregg14 + 100*dlambdachig14 - 27*lmMt) + 30*(40
      *dg1g1 - 40*dvg12 + 3*lmMt*(-3 + Xi))*pow2(c2beta) - 9*(-116 + 45*lmhtreeMt
      - 60*lmMt + 10*lmwMt + 15*sqrt(3.)*Pi + 20*asin(mhtree/(2.*p.MW))*sqrt(-1
      + (4*pow2(p.MW))/pow2(mhtree)))*pow4(c2beta) - 45*(-2 + lmzMt
      + 2*asin(std::abs(c2beta)/2.)*sqrt(-1 + 4/pow2(c2beta)))*(12 - 4
      *pow2(c2beta) + pow4(c2beta))))/4000.),
      "dmh2g14");
   const double dmh2g24 = isNaN(orderMap.at(EFTOrders::G24)*(-(v2*(-120 - 80
      *dlambdachig24 - 80*dlambdag24 - 80*dlambdaregg24 + 180*lmMt + 80*dg2g2
      *pow2(c2beta) - 80*dvg22*pow2(c2beta) - 90*lmMt*pow2(c2beta) + 30*lmMt*Xi
      *pow2(c2beta) + 6*pow4(c2beta) + 60*lmMt*pow4(c2beta) - 15*(-6 + 3
      *lmhtreeMt + sqrt(3)*Pi)*pow4(c2beta) - 5*(-2 + lmzMt + 2
      *asin(std::abs(c2beta)/2.)*sqrt(-1 + 4/pow2(c2beta)))*(12 - 4*pow2(c2beta)
      + pow4(c2beta)) - 10*(-2 + lmwMt + 2*asin(mhtree/(2.*p.MW))*sqrt(-1 + (4
      *pow2(p.MW))/pow2(mhtree)))*(12 - 4*pow2(c2beta) + pow4(c2beta))))/160.),
      "dmh2g24");
   const double dmh2g12yb2 = isNaN(orderMap.at(EFTOrders::G12YB2)*((v2*(10
      *dlambdayb2g12 + 3*(-6 + 2*dvyb2 + 3*lmhtreeMt - 3*lmMt)*pow2(c2beta))
      *pow2(cbeta))/20.),
      "dmh2g12yb2");
   const double dmh2g22yb2 = isNaN(orderMap.at(EFTOrders::G22YB2)*((v2*(2
      *dlambdayb2g22 + (-6 + 2*dvyb2 + 3*lmhtreeMt - 3*lmMt)*pow2(c2beta))
      *pow2(cbeta))/4.),
      "dmh2g22yb2");
   const double dmh2yb4 = isNaN(orderMap.at(EFTOrders::YB4)*(pow4(cbeta)*v2
      *(dlambdayb4 + 12*(2 - lmhtreeMt + lmMt))/2.),
      "dmh2yb4");
   const double dmh2g12yt2 = isNaN(orderMap.at(EFTOrders::G12YT2)*(pow2(sbeta)*v2*(10
      *dlambdayt2g12 + pow2(c2beta)*(6 + 6*dvyt2 - 9*lmMt))/20.),
      "dmh2g12yt2");
   const double dmh2g22yt2 = isNaN(orderMap.at(EFTOrders::G22YT2)*(pow2(sbeta)*v2*(2
      *dlambdayt2g22 + pow2(c2beta)*(2 + 2*dvyt2 - 3*lmMt))/4.),
      "dmh2g22yt2");
   const double dmh2yt4 = isNaN(orderMap.at(EFTOrders::YT4)*(pref_at*(12 * lmMt +
      thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT,
         RenSchemes::DRBARPRIME, omitMSSMLogs))),
      "dmh2yt4");
   const double dmh2g12ytau2 = isNaN(orderMap.at(EFTOrders::G12YTAU2)*(((10*dlambdaytau2g12
      + 3*(-2 + 2*dvytau2 + lmhtreeMt - lmMt)*pow2(c2beta))*pow2(cbeta)*v2)/20.),
      "dmh2g12ytau2");
   const double dmh2g22ytau2 = isNaN(orderMap.at(EFTOrders::G22YTAU2)*(((2*dlambdaytau2g22
      + (-2 + 2*dvytau2 + lmhtreeMt - lmMt)*pow2(c2beta))*pow2(cbeta)*v2)/4.),
      "dmh2g22ytau2");
   const double dmh2ytau4 = isNaN(orderMap.at(EFTOrders::YTAU4)*(pow4(cbeta)*v2*(8
      + dlambdaytau4 - 4*lmhtreeMt + 4*lmMt)/2.),
      "dmh2ytau4");

   // Loop factor
   const double k = 1/pow2(4.*Pi);

   return dmh2yt4 + k*(pow2(p.g1*p.g2)*dmh2g12g22 + pow4(p.g1)*dmh2g14 + pow4(p.g2)*
      dmh2g24 + pow2(p.g1*yb)*dmh2g12yb2 + pow2(p.g2*yb)*dmh2g22yb2 + pow4(yb)*
      dmh2yb4 + pow2(p.g1*ytau)*dmh2g12ytau2 + pow2(p.g2*ytau)*dmh2g22ytau2 +
      pow4(ytau)*dmh2ytau4 + pow2(p.g1*yt)*dmh2g12yt2 + pow2(p.g2*yt)*dmh2g22yt2);
}

/**
 *         Returns the 2-loop EFT contribution to the Higgs mass
 *         @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 *         @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 */
double Mh2EFTCalculator::getDeltaMh2EFT2Loop(int omitSMLogs, int omitMSSMLogs) const
{
   ThresholdCalculator thresholdCalculator(p, msq2);

   using std::log;
   const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
   // couplings
   const double v2 = pow2(p.vu) + pow2(p.vd);
   const double gt = sqrt(2)*p.Mt/std::sqrt(v2);
   const double g32 = pow2(p.g3);
   const double yt2 = pow2(sqrt(2.)*p.Mt/p.vu);
   const double yb2 = pow2(sqrt(2.)*p.Mb/p.vd);
   const double ytau2 = pow2(sqrt(2.)*p.Mtau/p.vd);
   const double yt4 = pow2(yt2);
   const double yb4 = pow2(yb2);
   const double yt6 = pow3(yt2);
   const double yb6 = pow3(yb2);
   const double ytau4 = pow2(ytau2);
   const double ytau6 = pow3(ytau2);
   const double beta = atan(p.vu/p.vd);
   const double tbeta = p.vu/p.vd;
   const double cbeta = cos(beta);
   const double c2beta = cos(2*beta);
   const double sbeta = sin(beta);
   const double Mhtree = std::abs(c2beta*p.MZ);
   const double lmhtreeMt = log(pow2(Mhtree / p.Mt));
   const double lmbMt = log(pow2(p.Mb / p.Mt));
   const double lmtauMt = log(pow2(p.Mtau / p.Mt));
   const double YBC = 1.;
   const double Xb = p.Ad(2,2) - p.mu*p.vu/p.vd;
   const double mQ32 = p.mq2(2,2);
   const double mD32 = p.md2(2,2);
   const double mD3 = sqrt(mD32);
   const double mQ3 = sqrt(mQ32);
   const double Xt = p.Au(2,2) - p.mu*p.vd/p.vu;
   const double gy = sqrt(3. / 5.)*p.g1;
   const double mU3 = sqrt(p.mu2(2,2));
   const double MR2 = pow2(p.scale);
   const double lmAMR = omitMSSMLogs*log(pow2(p.MA) / MR2);
   const double lmUMR = omitMSSMLogs*log(pow2(p.mu) / MR2);
   const double lm3MR = omitMSSMLogs*log(pow2(p.MG) / MR2);

   // 2-Loop prefactor at*as
   const double pref = 1./pow4(4*Pi) * pow2(p.Mt * gt * p.g3);

   // Threshold corrections
   const double dytas = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb4g32 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YB4_G32, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb4 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YB4, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb6 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YB6, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb2g12 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YB2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb2g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YB2_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayt4 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_AT, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dytyt = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_YT, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayt6 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YT6, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dvyt2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::VEV_YT2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaytau2g12 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaytau2g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dytauytau = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YTAU_YTAU, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaytau4 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU4, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdaytau6 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU6, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayt2yb4 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YT2_YB4, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayt2g12 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YT2_G12, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayt2g22 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YT2_G22, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayt4yb2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YT4_YB2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dytyb = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_YB, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dvytau2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::VEV_YTAU2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dvyb2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::VEV_YB2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dytauyb = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YTAU_YB, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb4ytau2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU2_YB4, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dlambdayb2ytau4 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::LAMBDA_YTAU4_YB2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dybyt = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YB_YT, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dybas = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YB_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);
   const double dybyb = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YB_YB, RenSchemes::DRBARPRIME, omitMSSMLogs);


   // Corrections to Mh
   const double dmh2yt4g32 = isNaN(orderMap.at(EFTOrders::G32YT4)*(pref*(96 * pow2(lmMt)
      + (-32 + 48 * dytas) * lmMt - 24 * dytas
      + thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT_AS,
         RenSchemes::DRBARPRIME, omitMSSMLogs))),
      "dmh2yt4g32");
   double dmh2yb4g32 = isNaN(orderMap.at(EFTOrders::G32YB4)*(((dlambdayb4g32 
      - 48*lmhtreeMt*(dybas + 16*lmMt) + 16*(3*dybas*(2 + lmMt)
      + 2*lmMt*(47 + 12*lmMt)))*v2*pow4(cbeta))/2.),
      "dmh2yb4g32");
   double dmh2yb6 = isNaN(orderMap.at(EFTOrders::YB6)*((v2*pow4(cbeta)*(
      -288*dybyb*(-2 + lmhtreeMt - lmMt)*pow2(Mhtree) + pow2(cbeta)*(-12*v2
      *(3*dlambdayb2g22*pow2(p.g2) + 5*dlambdayb2g12*pow2(gy)) + pow2(Mhtree)
      *(49 + 6*dlambdayb6 + 144*dvyb2 + 60*lmbMt - 144*dvyb2*lmhtreeMt - 5652*lmMt
      + 144*dvyb2*lmMt + 2592*lmhtreeMt*lmMt - 36*pow2(lmbMt) - 1296*pow2(lmMt) 
      + 6*pow2(Pi) + 6*dlambdayb4*(3 + (-6 + 2*dvyb2+ 3*lmhtreeMt - 3*lmMt)
      *pow4(YBC))))))/(12.*pow2(Mhtree))),
      "dmh2yb6");
   const double dmh2yt6= isNaN(orderMap.at(EFTOrders::YT6)*((v2*(4*dytyt*(
      -6 + dlambdayt4 + 12*lmMt) + (-12 + dlambdayt6 + dlambdayt4*(2 + 2*dvyt2
      - 3*lmMt) + 24*dvyt2*(-1 + lmMt) - 18*lmMt*(1 + 3*lmMt) - 2*pow2(Pi))
      *pow2(sbeta))*pow4(sbeta))/2.),
      "dmh2yt6");
   const double dmh2yb4ytau2 = isNaN(orderMap.at(EFTOrders::YTAU2YB4)*((v2*(
      -2*v2*(3*dlambdaytau2g22*pow2(p.g2) + 5*dlambdaytau2g12*pow2(gy))
      + lmhtreeMt*pow2(Mhtree)*(-24*dvytau2 + 96*lmMt + dlambdayb4*pow4(YBC)) 
      + pow2(Mhtree)*(dlambdayb4 + dlambdayb4ytau2 - 48*lmMt*(4 + lmMt) 
      - dlambdayb4*(2 + lmMt)*pow4(YBC) + 2*dvytau2*(12 + 12*lmMt 
      + dlambdayb4*pow4(YBC))))*pow6(cbeta))/(2.*pow2(Mhtree))),
      "dmh2yb4ytau2");
   const double dmh2yb2ytau4 = isNaN(orderMap.at(EFTOrders::YTAU4YB2)*((v2*(-2*v2
      *pow2(cbeta)*(3*dlambdayb2g22*pow2(p.g2) + 5*dlambdayb2g12*pow2(gy)) + 3
      *lmhtreeMt*(-16*dytauyb + (3*dlambdaytau4 - 8*dvyb2 + 24*lmMt)*pow2(cbeta))
      *pow2(Mhtree) + 3*(4*dytauyb*(8 + dlambdaytau4 + 4*lmMt) + (dlambdayb2ytau4 
      + 8*dvyb2*(1 + lmMt) - 6*lmMt*(9 + 2*lmMt) + dlambdaytau4*(2*dvyb2 - 3*(1 
      + lmMt)))*pow2(cbeta))*pow2(Mhtree))*pow4(cbeta))/(6.*pow2(Mhtree))),
      "dmh2yb2ytau4");
   const double dmh2ytau6 = isNaN(orderMap.at(EFTOrders::YTAU6)*(-(v2*(-12
      *dytauytau*(8 + dlambdaytau4 - 4*lmhtreeMt + 4*lmMt)*pow2(Mhtree) 
      + pow2(cbeta)*(2*v2*(3*dlambdaytau2g22*pow2(p.g2) + 5*dlambdaytau2g12
      *pow2(gy)) + pow2(Mhtree)*(6 - 3*dlambdaytau6 + 180*lmMt - 60*lmhtreeMt
      *lmMt + dlambdaytau4*(3 - 3*lmhtreeMt + 3*lmMt) - 6*dvytau2*(4 
      + dlambdaytau4 - 4*lmhtreeMt + 4*lmMt) - 21*lmtauMt + 30*pow2(lmMt) 
      + 9*pow2(lmtauMt) + pow2(Pi))))*pow4(cbeta))/(6.*pow2(Mhtree))),
      "dmh2ytau6");
   double dmh2yt2yb4 = isNaN(orderMap.at(EFTOrders::YT2YB4)*((v2*pow4(cbeta)
      *(96*dybyt*(2 + lmMt)*pow2(Mhtree) - 24*lmhtreeMt*pow2(Mhtree)*(4*dybyt 
      + (2*dvyt2 - 9*lmMt)*pow2(sbeta)) + pow2(sbeta)*(-4*v2*(3*dlambdayt2g22
      *pow2(p.g2) + 5*dlambdayt2g12*pow2(gy)) + pow2(Mhtree)*(2*dlambdayt2yb4 
      + 3*(8*lmbMt + 16*dvyt2*(1 + lmMt) - 3*(5 + 4*lmMt*(9 + 4*lmMt)) 
      - 2*pow2(Pi)) + 2*dlambdayb4*(2 + 2*dvyt2 - 3*lmMt)*pow4(YBC)))))
      /(4.*pow2(Mhtree))),
      "dmh2yt2yb4");
   const double dmh2yt4yb2 = isNaN(orderMap.at(EFTOrders::YB2YT4)*((v2*(4
      *dytyb*(-6 + dlambdayt4 + 12*lmMt) + pow2(cbeta)*(dlambdayt4yb2 
      + dlambdayt4*(-3 + 2*dvyb2 + 3*lmhtreeMt - 3*lmMt) + 3*(8*dvyb2*(-1 
      + lmMt) - 26*lmMt + 12*lmhtreeMt*lmMt - 6*pow2(lmMt) + pow2(Pi))))
      *pow4(sbeta))/2.),
      "dmh2yt4yb2");

   // the second term comes from the tanb resummation
   dmh2yb6 += isNaN(orderMap.at(EFTOrders::YB6)*((dlambdayb4*(2*Xb
      *F5(mQ3/mD3)*pow2(cbeta) - mD3*mQ3*(6*lmUMR + 4*F6(mD3/p.mu) + 8*F6(mQ3/p.mu)
      - 3*pow2(sbeta) + 6*lmAMR*pow2(sbeta))))/(8.*mD3*mQ3)),
      "dmh2yb6");
   dmh2yt2yb4 += isNaN(orderMap.at(EFTOrders::YT2YB4)*((dlambdayb4
      *(4*p.mu*Xt*F5(mQ3/mU3)*pow2(sbeta) + mQ3*mU3*(-8*p.mu*F6(mU3/p.mu) + p.mu*(9
      - 7*c2beta - 4*lmUMR + 2*(-5 + 3*c2beta)*lmAMR) - 8*Xt
      *F9(mQ3/p.mu,mU3/p.mu)*tbeta)))/(16.*mQ3*p.mu*mU3)),
      "dmh2yt2yb4");
   dmh2yb4g32 += isNaN(orderMap.at(EFTOrders::G32YB4)*((-4
      *dlambdayb4*(p.MG + p.MG*F6(mD3/p.MG) + p.MG*F6(mQ3/p.MG) - Xb*F9(mQ3/p.MG,mD3/p.MG)
      + p.MG*lm3MR))/(3.*p.MG)),
      "dmh2yb4g32");

   // Loop factor
   const double k2 = 1/pow4(4.*Pi);

   return k2*(g32*yb4*dmh2yb4g32 + yb6*dmh2yb6 + yt6*dmh2yt6 + yb4*ytau2
      *dmh2yb4ytau2 + yb2*ytau4*dmh2yb2ytau4 + ytau6*dmh2ytau6 + yt2*yb4
      *dmh2yt2yb4 + yt4*yb2*dmh2yt4yb2) + dmh2yt4g32;
}

/**
 *         Returns the 3-loop EFT contribution to the Higgs mass
 *         @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 *         @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 *         @param omitDeltaLambda3L an integer flag to disable the MSSM contribution to delta_lambda_3L
 */
double Mh2EFTCalculator::getDeltaMh2EFT3Loop(int omitSMLogs,
                                             int omitMSSMLogs,
                                             int omitDeltaLambda3L) const
{
   ThresholdCalculator thresholdCalculator(p, msq2);

   using std::log;

   const double catas2 = 248.1215180432007;
   const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
   // threshold correction of yt_as DRbar'
   const double dytas = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);
   // threshold correction of yt_as2 DRbar'
   const double dytas2 = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_AS2, RenSchemes::DRBARPRIME, omitMSSMLogs);
   // threshold correction of g3_as DRbar'
   const double dg3as = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::G3_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);

   const double gt = sqrt(2)*p.Mt/std::sqrt(pow2(p.vu) + pow2(p.vd));

   // 3-Loop prefactor at*as^2
   const double pref = 1./pow6(4*Pi) * pow2(p.Mt * gt * pow2(p.g3));

   return pref*(736 * pow3(lmMt) + (160 + 192 * dg3as + 384 * dytas) * pow2(lmMt)
      + (-128 * zt3 - 2056 / 3. + -64 * dg3as - 512 * dytas + 72 * pow2(dytas)
      + 48 * dytas2) * lmMt + 64 * dytas - 84 * pow2(dytas) - 24 * dytas2
      + catas2
      + omitDeltaLambda3L*thresholdCalculator.getThresholdCorrection(
         ThresholdVariables::LAMBDA_AT_AS2,
         RenSchemes::DRBARPRIME, omitMSSMLogs));
}

/**
 *   Returns the matching relation of delta_Lambda 3L for the degenerate mass case
 *   @param scale the renormalization scale
 *   @param mst1 the mass of the light stop quark
 *   @return delta_Lambda 3L
 */
double Mh2EFTCalculator::getDeltaLambdaDegenerate(double scale, double mst1, double Xt, int omitlogs) const
{
   using std::log;

   const double LS = omitlogs*log(pow2(scale / p.MSt(0))) + log(pow2(p.MSt(0) / mst1));

   const double gt = sqrt(2)*p.Mt/std::sqrt(pow2(p.vu) + pow2(p.vd));

   // 3-Loop prefactor
   const double pref = 1./pow6(4*Pi) * pow2(p.Mt * gt * pow2(p.g3));

   // to obtain delta_lambda one has to divide the difference of the two calculations by v^2
   const double v2 = pow2(p.vd) + pow2(p.vu);

   const double xt = Xt/mst1;
   const double catas2 = 248.1215180432007;

   const double deltaLambda3L = 2/27.*pref*(6082 - 27832*LS + 14856*pow2(LS)
      - 4032*pow3(LS) - 15408*zt3 + 1728*zt3*LS - 27*catas2/2.
      + xt*(7616*LS - 11712*pow2(LS) + 32*(-940 + 477*zt3))
      + pow2(xt)*(28848 - 2640*LS + 1008*pow2(LS) - 11880*zt3)
      + pow3(xt)*(160*LS + 864*pow2(LS) + 8*(2722 - 2259*zt3)))/v2;

   return deltaLambda3L;
}

} // namespace mh2_eft
} // namespace himalaya

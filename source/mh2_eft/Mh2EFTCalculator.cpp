#include "Mh2EFTCalculator.hpp"
#include "ThresholdCalculator.hpp"
#include "Hierarchies.hpp"
#include "Logger.hpp"
#include <iostream>

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

} // anonymous namespace
} // namespace mh2_eft
} // namespace himalaya

/**
 *	Constructor
 * 	@param p_ a HimalayaInterface struct
 * 	@param msq2_ the averaged squark mass of the first two generations squared
 * 	@param verbose a bool enable the output of the parameter validation. Enabled by default
 */
himalaya::mh2_eft::Mh2EFTCalculator::Mh2EFTCalculator(
   const himalaya::Parameters& p_, double msq2_, bool verbose)
   : p(p_), msq2(msq2_)
{
   p.validate(verbose);;

   if (!std::isfinite(msq2_))
      msq2 = p.calculateMsq2();
   
   // fill order map
   orderMap.clear();
   for (int i = EFTOrders::FIRST; i < EFTOrders::NUMBER_OF_EFT_ORDERS; i++) {
      orderMap.emplace(i, 1u);
   }
}

void himalaya::mh2_eft::Mh2EFTCalculator::setCorrectionFlag(int variable, int enable){
   if(enable < 0 || enable > 1) INFO_MSG("You can only enable (1) or disable (0) corrections!");
   if(orderMap.find(variable) == orderMap.end()) INFO_MSG("Your variable is not defined in the EFTOrders enum!");
   orderMap.at(variable) = enable;
}


/**
 * 	Returns the 1-loop EFT contribution to the Higgs mass
 * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 */
double himalaya::mh2_eft::Mh2EFTCalculator::getDeltaMh2EFT1Loop(int omitSMLogs, 
								int omitMSSMLogs){
   ThresholdCalculator thresholdCalculator(p, msq2);
   
   using std::log;
   const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
   
   const double gt = sqrt(2)*p.Mt/std::sqrt(pow2(p.vu) + pow2(p.vd));
   
   // 1-Loop prefactor at
   const double pref_at = 1./pow2(4*Pi) * pow2(p.Mt * gt);
   
   const double v2 = pow2(p.vu) + pow2(p.vd);
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
   const int Xi = 1;	// gauge parameter
   
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

   // corrections to Mh2
   const double dmh2g12g22 = orderMap.at(EFTOrders::G12G22)*(pow2(cbeta)*v2*(10
      *dlambdayb2g12 - 9*pow2(c2beta)*(2 + lmMt - lmhtreeMt))/20.);
   const double dmh2g14 = orderMap.at(EFTOrders::G14)*((v2*(20*(-36 
      + 100*dlambdag14 + 100*dlambdaregg14 + 100*
      dlambdachig14 - 27*lmMt + 27*lmzMt) - 30*(40*dg1g1 + 6*(-2 + lmzMt) +
      3*lmMt*(-3 + Xi))*pow2(c2beta) + 9*(-126 + 45*lmhtreeMt - 60*lmMt + 10*
      lmwMt + 5*lmzMt + 15*sqrt(3)*Pi)*pow4(c2beta) + 180*asin(mhtree/(2*p.MW)
      )*sqrt(-1 + (4*pow2(p.MW))/pow2(mhtree))*pow4(c2beta) + 90*asin(c2beta/
      2.)*sqrt(-1 + 4/pow2(c2beta))*(12 - 4*pow2(c2beta) + pow4(c2beta))))/4000.);
   const double dmh2g24 = orderMap.at(EFTOrders::G24)*((v2*(20*(-12 
      + 4*dlambdachig24 + 4*dlambdag24 + 4*
      dlambdaregg24 - 9*lmMt + 6*lmwMt + 3*lmzMt) - 10*(8*dg2g2 + 2*(-6 + 
      2*lmwMt + lmzMt) + 3*lmMt*(-3 + Xi))*pow2(c2beta) + (-126 + 45*lmhtreeMt -
      60*lmMt + 10*lmwMt + 5*lmzMt + 15*sqrt(3.)*Pi)*pow4(c2beta) + 10*asin(
      c2beta/2.)*sqrt(-1 + 4/pow2(c2beta))*(12 - 4*pow2(c2beta) + pow4(c2beta))
      + 20*asin(mhtree/(2*p.MW))*sqrt(-1 + (4*pow2(p.MW))/pow2(mhtree))*(12 - 4*
      pow2(c2beta) + pow4(c2beta))))/160.);
   const double dmh2g12yb2 = orderMap.at(EFTOrders::G12YB2)*((v2*(20*(-12 
      + 10*dlambdag12g22 + 10*dlambdaregg12g22 + 10*
      dlambdachig12g22 - 9*lmMt + 9*lmzMt) - 60*(-4 + lmwMt + lmzMt + lmMt*(-
      3 + Xi))*pow2(c2beta) + 60*asin(mhtree/(2*p.MW))*(-2 + pow2(c2beta))*
      pow2(c2beta)*sqrt(-1 + (4*pow2(p.MW))/pow2(mhtree)) + 3*(-126 + 45*
      lmhtreeMt - 60*lmMt + 10*lmwMt + 5*lmzMt + 15*sqrt(3.)*Pi)*pow4(c2beta)
      + 30*asin(c2beta/2.)*sqrt(-1 + 4/pow2(c2beta))*(12 - 4*pow2(c2beta) +
      pow4(c2beta))))/400.);
   const double dmh2g22yb2 = orderMap.at(EFTOrders::G22YB2)*(pow2(cbeta)*v2*(2
      *dlambdayb2g22 + 3*pow2(c2beta)*(-2 + lmhtreeMt - lmMt))/4.);
   const double dmh2yb4 = orderMap.at(EFTOrders::YB4)*(pow4(cbeta)*v2
      *(dlambdayb4 + 12*(2 - lmhtreeMt + lmMt))/2.);
   const double dmh2g12yt2 = orderMap.at(EFTOrders::G12YT2)*(pow2(sbeta)*v2*(10
      *dlambdayt2g12 + pow2(c2beta)*(6 + 6*dvyt2 - 9*lmMt))/20.);
   const double dmh2g22yt2 = orderMap.at(EFTOrders::G22YT2)*(pow2(sbeta)*v2*(2
      *dlambdayt2g22 + pow2(c2beta)*(2 + 2*dvyt2 - 3*lmMt))/4.);
   const double dmh2yt4 = orderMap.at(EFTOrders::YT4)*(pref_at*(12 * lmMt +
      thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT,
	 RenSchemes::DRBARPRIME, omitMSSMLogs)));
   const double dmh2g12ytau2 = orderMap.at(EFTOrders::G12YTAU2)*(pow2(cbeta)*v2*(
      10*dlambdaytau2g12 + 3*pow2(c2beta)*(-2 + lmhtreeMt - lmMt))/20.);
   const double dmh2g22ytau2 = orderMap.at(EFTOrders::G22YTAU2)*(pow2(cbeta)*v2*(2
      *dlambdaytau2g22 + pow2(c2beta)*(-2 + lmhtreeMt - lmMt))/4.);
   const double dmh2ytau4 = orderMap.at(EFTOrders::YTAU4)*(pow4(cbeta)*v2*(8
      + dlambdaytau4 - 4*lmhtreeMt + 4*lmMt)/2.);

   // Loop factor
   const double k = 1/pow2(4.*Pi);

   return dmh2yt4 + k*(pow2(p.g1*p.g2)*dmh2g12g22 + pow4(p.g1)*dmh2g14 + pow4(p.g2)*
      dmh2g24 + pow2(p.g1*yb)*dmh2g12yb2 + pow2(p.g2*yb)*dmh2g22yb2 + pow4(yb)*
      dmh2yb4 + pow2(p.g1*ytau)*dmh2g12ytau2 + pow2(p.g2*ytau)*dmh2g22ytau2 + 
      pow4(ytau)*dmh2ytau4 + pow2(p.g1*yt)*dmh2g12yt2 + pow2(p.g2*yt)*dmh2g22yt2);
}

/**
 * 	Returns the 2-loop EFT contribution to the Higgs mass
 * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 */
double himalaya::mh2_eft::Mh2EFTCalculator::getDeltaMh2EFT2Loop(int omitSMLogs,
								int omitMSSMLogs){
   ThresholdCalculator thresholdCalculator(p, msq2);
   
   using std::log;
   const double lmMt = omitSMLogs * log(pow2(p.scale / p.Mt));
   // couplings
   const double gt = sqrt(2)*p.Mt/std::sqrt(pow2(p.vu) + pow2(p.vd));
   /*const double g32 = pow2(p.g3);
   const double yt2 = pow2(sqrt(2.)*p.Mt/p.vu);
   const double yb2 = pow2(sqrt(2.)*p.Mb/p.vd);
   const double ytau2 = pow2(sqrt(2.)*p.Mtau/p.vd);
   const double yt4 = pow2(yt2);
   const double yb4 = pow2(yb2);
   const double yt6 = pow3(yt2);
   const double yb6 = pow3(yb2);
   const double ytau6 = pow3(ytau2);
   const double v2 = pow2(p.vu) + pow2(p.vd);
   const double beta = atan(p.vu/p.vd);
   const double cbeta = cos(beta);
   const double c2beta = cos(2*beta);
   const double sbeta = sin(beta);
   const double mhtree = std::abs(c2beta*p.MZ);
   const double lmhtreeMt = log(pow2(mhtree / p.Mt));
   const double lmbMt = log(pow2(p.Mb / p.Mt));
   const double lmtauMt = log(pow2(p.Mtau / p.Mt));*/
   
   // 2-Loop prefactor at*as
   const double pref = 1./pow4(4*Pi) * pow2(p.Mt * gt * p.g3);
   
   // Threshold corrections
   const double dytas = thresholdCalculator.getThresholdCorrection(
      ThresholdVariables::YT_AS, RenSchemes::DRBARPRIME, omitMSSMLogs);
   /*const double dlambdayb4g32 = thresholdCalculator.getThresholdCorrection(
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
      ThresholdVariables::YT_YB, RenSchemes::DRBARPRIME, omitMSSMLogs);*/
   
   // Corrections to Mh
   const double dmh2yt4g32 = orderMap.at(EFTOrders::G32YT4)*(pref*(96 * pow2(lmMt)
      + (-32 + 48 * dytas) * lmMt - 24 * dytas
      + thresholdCalculator.getThresholdCorrection(ThresholdVariables::LAMBDA_AT_AS,
	 RenSchemes::DRBARPRIME, omitMSSMLogs)));
   /*const double dmh2yb4g32 = orderMap.at(EFTOrders::G32YB4)*(pow4(cbeta)
      *(dlambdayb4g32 + 32*lmMt*(47 - 24*lmhtreeMt + 12*lmMt))*v2/2.);
   const double dmh2yb6 = orderMap.at(EFTOrders::YB6)*(pow6(cbeta)*v2*(29 + 2
      *dlambdayb6 + 4*lmbMt*(-38 + 9*lmbMt) + 6*dlambdayb4*(-1 + lmhtreeMt - lmMt)
      + 12*(-257 + 72*lmhtreeMt - 36*lmMt)*lmMt - 6*pow2(Pi) - 12*(dlambdayb2g12
      *pow2(p.g1) + dlambdayb2g22*pow2(p.g2))*v2/pow2(mhtree))/4.);
   const double dmh2yt6= orderMap.at(EFTOrders::YT6)*(pow4(sbeta)*(24*dytyt*(-6
      + dlambdayt4 + 12*lmMt) + 113*pow2(sbeta) + 6*(dlambdayt6 + dlambdayt4*(2
      + 2*dvyt2 - 3*lmMt) - 3*(-2*pow2(lmbMt) - 8*dvyt2*(lmMt - 1) + 6*lmMt*(1
      + 3*lmMt) + pow2(Pi))*pow2(sbeta))*v2)/12.);
   const double dmh2yb4ytau2 = orderMap.at(EFTOrders::YTAU2YB4)*(pow6(cbeta)*v2
      *(dlambdayb4*(lmhtreeMt - 1 - lmMt)*pow2(mhtree) - 48*lmMt*(4 - 2*lmhtreeMt
      + lmMt)*pow2(mhtree) - 6*(dlambdaytau2g12*pow2(p.g1) + dlambdaytau2g22
      *pow2(p.g2)))*v2/(2*pow2(mhtree)));
   const double dmh2yt4ytau2 = 0.;	//??
   const double dmh2ytau6 = orderMap.at(EFTOrders::YTAU6)*(pow4(cbeta)*v2*(12
      *dytauytau*(8 + dlambdaytau4 - 4*lmhtreeMt + 4*lmMt) + pow2(cbeta)*(3
      *dlambdaytau6 + 3*dlambdaytau4*(lmhtreeMt - 1 - lmMt) - 2*(6 + 15*lmMt
      *(6 - 2*lmhtreeMt + lmMt) + 3*lmtauMt*(3*lmtauMt - 7) + pow2(Pi)) - 6
      *(dlambdaytau2g12*pow2(p.g1) + dlambdaytau2g22*pow2(p.g2))*v2
      /pow2(mhtree)))/6.);
   const double dmh2yt2yb4 = orderMap.at(EFTOrders::YT2YB4)*(pow4(cbeta)*pow2(
      sbeta)*v2*(pow2(mhtree)*(dlambdayt2yb4 + dlambdayb4*(2 + 2*dvyt2 - 3*lmMt)
      + 6*(-2*lmbMt + 3*(-9 + 6*lmhtreeMt - 4*lmMt)*lmMt + 4*dvyt2*(1 - lmMt
      - lmMt) + 2*(5 + pow2(Pi)))) - 6*(dlambdayt2g12*pow2(p.g1) + dlambdayt2g22
      *pow2(p.g2)*v2))/(2.*pow2(mhtree)));
   const double dmh2yt4yb2 = orderMap.at(EFTOrders::YB2YT4)*(pow4(sbeta)*v2*(4
      *dytyb*(-6 + dlambdayt4 + 12*lmMt) + pow2(cbeta)*(dlambdayt4yb2 + 3
      *dlambdayt4*(lmhtreeMt - 1 - lmMt) - 6*(5 + lmMt*(13 - 6*lmhtreeMt + 3
      *lmMt) + 2*pow2(Pi))))/2.);*/
   
   // Loop factor
   //const double k2 = 1/pow4(4.*Pi);

   return /*k2*(g32*yb4*dmh2yb4g32 + yb6*dmh2yb6 + yt6*dmh2yt6 + yb4*ytau2
      *dmh2yb4ytau2 + yt4*ytau2*dmh2yt4ytau2 + ytau6*dmh2ytau6 + yt2*yb4
      *dmh2yt2yb4 + yt4*yb2*dmh2yt4yb2)*/ + dmh2yt4g32;
}

/**
 * 	Returns the 3-loop EFT contribution to the Higgs mass
 * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
 * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
 * 	@param omitDeltaLambda3L an integer flag to disable the MSSM contribution to delta_lambda_3L
 */
double himalaya::mh2_eft::Mh2EFTCalculator::getDeltaMh2EFT3Loop(int omitSMLogs,
								int omitMSSMLogs,
								int omitDeltaLambda3L){
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
double himalaya::mh2_eft::Mh2EFTCalculator::getDeltaLambdaDegenerate(double scale, double mst1, double Xt, int omitlogs) const{
   using std::log;
   
   const double li4 = 0.5174790616738993; // Li4(1/2)
   
   const double LS = omitlogs*log(pow2(scale / p.MSt(0))) + log(pow2(p.MSt(0) / mst1));
   
   const double gt = sqrt(2)*p.Mt/std::sqrt(pow2(p.vu) + pow2(p.vd));
   
   // 3-Loop prefactor
   const double pref = 1./pow6(4*Pi) * pow2(p.Mt * gt * pow2(p.g3));
   
   // to obtain delta_lambda one has to divide the difference of the two calculations by v^2
   const double v2 = pow2(p.vd) + pow2(p.vu);
   
   // the first 16 is a relict due to the zeta parameterization in the first draft
   const double deltaLambda3L = 16.*pref*((29365 + 23040*li4 - 49320*zt3 + 160*LS*
	(-226 + 27*zt3) + 26520*pow2(LS) - 80*Xt*(823 - 100*LS - 477*zt3 + 366*
	pow2(LS)) + 30*(-67 - 220*LS - 990*zt3 + 84*pow2(LS))*pow2(Xt) - 960*
	pow2(Pi)*pow2(log2) - 10080*pow3(LS)+ 20*(3568 + 20*LS - 2259*zt3 + 108*
	pow2(LS))*pow3(Xt) - 176*pow4(Pi) + 960*pow4(log2))/540.)/v2;
   return deltaLambda3L;
}

/**
 * 	A function which maps a boolean to a string.
 * 	@param tf a boolean.
 * 	@return A string which is 'true' if tf is true or 'false' if tf is false.
 */
std::string himalaya::mh2_eft::Mh2EFTCalculator::tf(const bool tf){
   return tf ? "\033[1;32mtrue\033[0m" : "\033[1;31mfalse\033[0m";
}

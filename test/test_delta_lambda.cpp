#include "doctest.h"
#include "HierarchyCalculator.hpp"
#include "Mh2EFTCalculator.hpp"


#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

namespace {

double pow2(double x) { return x*x; }

himalaya::Parameters make_point()
{
   const double MS = 100000.;
   const double xt = 2.;
   const double tb = 5.;
   const double beta = std::atan(tb);
   const double v = 245.;

   const double msq2 = MS;
   const double msu2 = MS;
   const double mg   = MS;
   const double msl2 = MS;
   const double mse2 = MS;

   const double Xt = xt * std::sqrt(std::sqrt(msq2 * msu2));

   himalaya::Parameters pars;
   pars.scale = MS;
   pars.mu    = MS;
   pars.g1    = 0.448494;
   pars.g2    = 0.607902;
   pars.g3    = 1.05733;
   pars.vu    = v*std::sin(beta);
   pars.vd    = v*std::cos(beta);
   pars.mq2   << pow2(msq2), 0, 0, 0, pow2(msq2), 0, 0, 0, pow2(msq2);
   pars.md2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS);
   pars.mu2   << pow2(msu2), 0, 0, 0, pow2(msu2), 0, 0, 0, pow2(msu2);
   pars.ml2   << pow2(msl2), 0, 0, 0, pow2(msl2), 0, 0, 0, pow2(msl2);
   pars.me2   << pow2(mse2), 0, 0, 0, pow2(mse2), 0, 0, 0, pow2(mse2);
   pars.Au    << 0, 0, 0, 0, 0, 0, 0, 0, Xt + pars.mu/tb;
   pars.Ae    << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Ad    << 0, 0, 0, 0 ,0 ,0 ,0 ,0 ,0;
   pars.M1    = MS;
   pars.M2    = MS;
   pars.MG    = mg;
   pars.MW    = 74.597;
   pars.MZ    = 85.7704;
   pars.Mt    = 144.337;
   pars.Mb    = 2.37054;
   pars.Mtau  = 1.2;
   pars.MA    = MS;

   pars.validate(false);

   return pars;
}

// point with small Xt, because of missing higher order Xt terms in H3m
himalaya::Parameters make_degenerate_point(double eps = 0.)
{
   // MS = 1 TeV, Xt = 0.1 MS, TB = 5
   himalaya::Parameters pars;
   pars.scale = 1000;
   pars.mu    = 1000;
   pars.g3    = 1.05733;
   pars.vd    = 45.8309;
   pars.vu    = 229.155;
   pars.mq2   << 1e+06, 0, 0, 0, 1e+06, 0, 0, 0, 1e+06*(1 + eps);
   pars.md2   << 1e+06, 0, 0, 0, 1e+06, 0, 0, 0, 1e+06;
   pars.mu2   << 1e+06, 0, 0, 0, 1e+06, 0, 0, 0, 1e+06*(1 - eps);
   pars.Au    << 0, 0, 0, 0, 0, 0, 0, 0, 300;
   pars.Ad    << 0, 0, 0, 0, 0, 0, 0, 0, 5000;
   pars.MG    = 1000*(1 - eps);
   pars.MW    = 74.597;
   pars.MZ    = 85.7704;
   pars.Mt    = 144.337;
   pars.Mb    = 2.37054;
   pars.MA    = 1000;

   return pars;
}

/// calculates lightest mass eigenvalue of given mass matrix
double calc_Mh2(const Eigen::Matrix2d& Mh_mat)
{
   const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(Mh_mat);
   return es.eigenvalues()(0);
}

/// calculates Mh^2 in the EFT at tree level
double calc_Mh2_EFT_0L(const himalaya::Parameters& pars)
{
   const double beta = std::atan(pars.vu/pars.vd);
   const double cos_2beta = std::cos(2*beta);

   return pow2(pars.MZ * cos_2beta);
}

/// calculates Mh^2 in the EFT at 1-loop level
double calc_Mh2_EFT_1L(const himalaya::Parameters& pars)
{
   himalaya::mh2_eft::Mh2EFTCalculator mhc(pars);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G12G22, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G12YB2, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G14, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G24, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G12YB2, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G22YB2, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::YB4, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G12YTAU2, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G22YTAU2, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::YTAU4, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G12YT2, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G22YT2, 0);

   return mhc.getDeltaMh2EFT1Loop(1,1);
}

/// calculates Mh^2 in the EFT at 2-loop level
double calc_Mh2_EFT_2L(const himalaya::Parameters& pars)
{
   himalaya::mh2_eft::Mh2EFTCalculator mhc(pars);
   mhc.setCorrectionFlag(himalaya::EFTOrders::G32YB4, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::YB6, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::YT6, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::YTAU2YB4, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::YTAU6, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::YT2YB4, 0);
   mhc.setCorrectionFlag(himalaya::EFTOrders::YB2YT4, 0);

   return mhc.getDeltaMh2EFT2Loop(1,1);
}

/// calculates Mh^2 in the EFT at 3-loop level
/// zeta_switch == 0: all logs (SM + MSSM) + zeta w/o any logs
/// zeta_switch == 1: SM logs + zeta w/ MSSM logs from H3m
/// zeta_switch == 2: SM logs + zeta w/ MSSM logs from EFT
double calc_Mh2_EFT_3L(const himalaya::Parameters& pars, int zeta_switch)
{
   auto hc = himalaya::HierarchyCalculator(pars);
   const auto ho = hc.calculateDMh3L(false);

   const double zeta_3L_const    = ho.getDLambdaNonLog();
   const double zeta_3L_eft      = ho.getDLambdaEFT();
   const double zeta_3L_himalaya = ho.getDLambdaH3m();
   const double v2 = pow2(pars.vu) + pow2(pars.vd);

   himalaya::mh2_eft::Mh2EFTCalculator mhc(pars);

   // calculate only logs
   const double DMh2_EFT_3L_logs_SM =
      mhc.getDeltaMh2EFT3Loop(1,0,4) - mhc.getDeltaMh2EFT3Loop(0,0,4);
   const double DMh2_EFT_3L_logs_all =
      mhc.getDeltaMh2EFT3Loop(1,1,4) - mhc.getDeltaMh2EFT3Loop(0,0,4);
   const double DMh2_EFT_3L_logs_mixed =
      + mhc.getDeltaMh2EFT3Loop(1,1,4) // all logs + const
      - mhc.getDeltaMh2EFT3Loop(1,0,4) // subtract SM logs + const
      - mhc.getDeltaMh2EFT3Loop(0,1,4) // subtract MSSM logs + const
      + mhc.getDeltaMh2EFT3Loop(0,0,4); // add constant

   double DMh2_3L = 0.;

   switch (zeta_switch) {
   case 0: DMh2_3L = DMh2_EFT_3L_logs_all + zeta_3L_const*v2;
      break;
   case 1: DMh2_3L = DMh2_EFT_3L_logs_SM + DMh2_EFT_3L_logs_mixed + zeta_3L_himalaya*v2;
      break;
   case 2: DMh2_3L = DMh2_EFT_3L_logs_SM + DMh2_EFT_3L_logs_mixed + zeta_3L_eft*v2;
      break;
   default: INFO("Error: unknow switch value : " << zeta_switch);
   }

   return DMh2_3L;
}

} // anonymous namespace

TEST_CASE("test_delta_lambda")
{
   using namespace himalaya;

   const auto pars = make_point();
   auto hc = HierarchyCalculator(pars);
   const auto ho = hc.calculateDMh3L(false);

   const auto DMh_0L = ho.getDMh(0);
   const auto DMh_1L = ho.getDMh(1);
   const auto DMh_2L = ho.getDMh(2);
   const auto DMh_3L = ho.getDMh(3);

   const auto Mh2_0L = calc_Mh2(DMh_0L);
   const auto Mh2_1L = calc_Mh2(DMh_0L + DMh_1L);
   const auto Mh2_2L = calc_Mh2(DMh_0L + DMh_1L + DMh_2L);
   const auto Mh2_3L = calc_Mh2(DMh_0L + DMh_1L + DMh_2L + DMh_3L);

   // 3L expansion uncertainty
   const double Mh2_3L_uncert = ho.getDMhExpUncertainty(3);
   const double Mh2_3L_uncert_rel = Mh2_3L_uncert / Mh2_3L;

   const double Mh_0L = std::sqrt(Mh2_0L);
   const double Mh_1L = std::sqrt(Mh2_1L);
   const double Mh_2L = std::sqrt(Mh2_2L);
   const double Mh_3L = std::sqrt(Mh2_3L);
   const double Mh_3L_uncert_percent = Mh2_3L_uncert_rel*100;

   INFO("Mh(0L) = " << Mh_0L << " GeV");
   INFO("Mh(1L) = " << Mh_1L << " GeV");
   INFO("Mh(2L) = " << Mh_2L << " GeV");
   INFO("Mh(3L) = (" << Mh_3L << " +- " << Mh2_3L_uncert << ") GeV"
        << " (" << Mh_3L_uncert_percent << "%)");

   CHECK_CLOSE(Mh2_0L, calc_Mh2_EFT_0L(pars), 1e-5);

   CHECK_CLOSE(Mh2_1L, calc_Mh2_EFT_0L(pars)
                     + calc_Mh2_EFT_1L(pars), 1e-6);

   CHECK_CLOSE(Mh2_2L, calc_Mh2_EFT_0L(pars)
                     + calc_Mh2_EFT_1L(pars)
                     + calc_Mh2_EFT_2L(pars), 1e-6);

   CHECK_CLOSE(Mh2_3L, calc_Mh2_EFT_0L(pars)
                     + calc_Mh2_EFT_1L(pars)
                     + calc_Mh2_EFT_2L(pars)
                     + calc_Mh2_EFT_3L(pars, 0), 1e-4);

   CHECK_CLOSE(Mh2_3L, calc_Mh2_EFT_0L(pars)
                     + calc_Mh2_EFT_1L(pars)
                     + calc_Mh2_EFT_2L(pars)
                     + calc_Mh2_EFT_3L(pars, 1), 1e-4);

   CHECK_CLOSE(Mh2_3L, calc_Mh2_EFT_0L(pars)
                     + calc_Mh2_EFT_1L(pars)
                     + calc_Mh2_EFT_2L(pars)
                     + calc_Mh2_EFT_3L(pars, 2), 1e-4);
}

TEST_CASE("test_lambda_limit_degenerate")
{
   using namespace himalaya;

   const auto pars = make_degenerate_point(0.03);
   auto hc = HierarchyCalculator(pars);
   const auto ho = hc.calculateDMh3L(false);

   const auto z2_gen_Himalaya = ho.getDLambdaH3m();
   const auto z2_gen_EFT      = ho.getDLambdaEFT();
   const auto uncertainty     = ho.getDLambdaUncertainty(3);

   INFO("uncertainty: " << uncertainty);

   CHECK_CLOSE(z2_gen_Himalaya, z2_gen_EFT, uncertainty);
}

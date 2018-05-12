#include "doctest.h"
#include "HierarchyCalculator.hpp"
#include "Mh2EFTCalculator.hpp"


#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

namespace {

const double Pi = 3.141592653589793;

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

   const double Xt = xt * std::sqrt(std::sqrt(msq2 * msu2));

   himalaya::Parameters pars;
   pars.scale = MS;
   pars.mu    = MS;
   pars.g3    = 1.05733;
   pars.vu    = v*std::sin(beta);
   pars.vd    = v*std::cos(beta);
   pars.mq2   << pow2(msq2), 0, 0, 0, pow2(msq2), 0, 0, 0, pow2(msq2);
   pars.md2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS);
   pars.mu2   << pow2(msu2), 0, 0, 0, pow2(msu2), 0, 0, 0, pow2(msu2);
   pars.At    = Xt + pars.mu/tb;
   pars.Ab    = 0;
   pars.MG    = mg;
   pars.MW    = 74.597;
   pars.MZ    = 85.7704;
   pars.Mt    = 144.337;
   pars.Mb    = 2.37054;
   pars.MA    = MS;

   pars.validate(false);

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

   return mhc.getDeltaMh2EFT1Loop(1,1);
}

/// calculates Mh^2 in the EFT at 2-loop level
double calc_Mh2_EFT_2L(const himalaya::Parameters& pars)
{
   himalaya::mh2_eft::Mh2EFTCalculator mhc(pars);

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

   const double zeta_3L_const    = ho.getDeltaLambdaNonLog();
   const double zeta_3L_eft      = ho.getDeltaLambdaEFT();
   const double zeta_3L_himalaya = ho.getDeltaLambdaHimalaya();

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
   case 0: DMh2_3L = DMh2_EFT_3L_logs_all + zeta_3L_const;
      break;
   case 1: DMh2_3L = DMh2_EFT_3L_logs_SM + DMh2_EFT_3L_logs_mixed + zeta_3L_himalaya;
      break;
   case 2: DMh2_3L = DMh2_EFT_3L_logs_SM + DMh2_EFT_3L_logs_mixed + zeta_3L_eft;
      break;
   default: INFO("Error: unknow switch value : " << zeta_switch);
   }

   return DMh2_3L;
}

} // anonymous namespace

TEST_CASE("test_lambda_normalization")
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
   const double Mh2_3L_uncert = ho.getExpUncertainty(3);
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

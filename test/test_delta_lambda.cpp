#include "doctest.h"

#include "himalaya/HierarchyCalculator.hpp"

#include "Mh2EFTCalculator.hpp"
#include "Flags.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

namespace {

double pow2(double x) { return x*x; }

himalaya::Parameters make_point(double MS = 100000., double tb = 5., double xt = 2.)
{
   const double beta = std::atan(tb);
   const double v = 245.;

   const double msq = MS;
   const double msu = MS;
   const double mg  = MS;
   const double msl = MS;
   const double mse = MS;

   const double Xt = xt * std::sqrt(msq * msu);

   himalaya::Parameters pars;
   pars.scale = MS;
   pars.mu    = MS;
   pars.g1    = 0.448494;
   pars.g2    = 0.607902;
   pars.g3    = 1.05733;
   pars.vu    = v*std::sin(beta);
   pars.vd    = v*std::cos(beta);
   pars.mq2   << pow2(msq), 0, 0, 0, pow2(msq), 0, 0, 0, pow2(msq);
   pars.md2   << pow2(MS) , 0, 0, 0, pow2(MS) , 0, 0, 0, pow2(MS);
   pars.mu2   << pow2(msu), 0, 0, 0, pow2(msu), 0, 0, 0, pow2(msu);
   pars.ml2   << pow2(msl), 0, 0, 0, pow2(msl), 0, 0, 0, pow2(msl);
   pars.me2   << pow2(mse), 0, 0, 0, pow2(mse), 0, 0, 0, pow2(mse);
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
   using namespace himalaya::mh2_eft;

   himalaya::mh2_eft::Mh2EFTCalculator mhc(pars);
   mhc.setCorrectionFlag(CouplingOrders::G12G22, 0);
   mhc.setCorrectionFlag(CouplingOrders::G12YB2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G14, 0);
   mhc.setCorrectionFlag(CouplingOrders::G24, 0);
   mhc.setCorrectionFlag(CouplingOrders::G12YB2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G22YB2, 0);
   mhc.setCorrectionFlag(CouplingOrders::YB4, 0);
   mhc.setCorrectionFlag(CouplingOrders::G12YTAU2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G22YTAU2, 0);
   mhc.setCorrectionFlag(CouplingOrders::YTAU4, 0);
   mhc.setCorrectionFlag(CouplingOrders::G12YT2, 0);
   mhc.setCorrectionFlag(CouplingOrders::G22YT2, 0);

   return mhc.getDeltaMh2EFT1Loop(1,1);
}

/// calculates Mh^2 in the EFT at 2-loop level
double calc_Mh2_EFT_2L(const himalaya::Parameters& pars)
{
   using namespace himalaya::mh2_eft;

   himalaya::mh2_eft::Mh2EFTCalculator mhc(pars);
   mhc.setCorrectionFlag(CouplingOrders::G32YB4, 0);
   mhc.setCorrectionFlag(CouplingOrders::YB6, 0);
   mhc.setCorrectionFlag(CouplingOrders::YT6, 0);
   mhc.setCorrectionFlag(CouplingOrders::YTAU2YB4, 0);
   mhc.setCorrectionFlag(CouplingOrders::YTAU6, 0);
   mhc.setCorrectionFlag(CouplingOrders::YT2YB4, 0);
   mhc.setCorrectionFlag(CouplingOrders::YB2YT4, 0);

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
                     + calc_Mh2_EFT_2L(pars), 1e-5);

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

/// Compare analytic expressions Eqs.(43)-(44) from arXiv:1807.03509
/// with Himalaya calculation.
TEST_CASE("test_lambda_degenerate")
{
   using namespace himalaya;

   const auto MS = 5000.;
   const auto tb = 20.;
   const auto xt = -2.;
   const auto pars = make_point(MS, tb, xt);

   const double PI = 3.141592653589793;
   const double z3 = 1.202056903159594; // Zeta[3]
   const double k = 0.006332573977646111; // 1/(4Pi)^2
   const double sqrt2 = 1.414213562373095; // Sqrt[2]
   const double Q = pars.scale;
   const double Q2 = pow2(Q);
   const double MS2 = pow2(MS);
   const double LS = std::log(Q2/MS2);
   const double LS2 = pow2(LS);
   const double LS3 = LS2*LS;
   const double xt2 = pow2(xt);
   const double xt3 = xt2*xt;
   const double yt = sqrt2*pars.Mt/pars.vu;
   const double g3 = pars.g3;
   const double alt = pow2(yt)/(4.*PI); // Eq.(7) from arXiv:1807.03509
   const double alt2 = pow2(alt);
   const double as = pow2(g3)*k; // Eq.(7) from arXiv:1807.03509
   const double as2 = pow2(as);
   const double sb = std::sin(std::atan(tb)); // sin(beta)
   const double sb4 = pow2(pow2(sb)); // sin^4(beta)

   auto hc = HierarchyCalculator(pars);
   const auto ho = hc.calculateDMh3L(false);

   // Delta lambda calculated with Himalaya
   const auto dl3 = ho.getDLambda(3);
   // delta lambda calculated with Himalaya (shift DR' -> MS)
   const auto dl3_shift = ho.getDLambdaDRbarPrimeToMSbarShift(3);

   /*
     Eq.(11) from arXiv:1807.03509

     cSM20 = (-1888/9 + 160 Zeta[3] + 7424/45 Zeta[2]^2 - 1024/3 PolyLog[4,1/2] - 512/9 PolyLog[2,1/2]^2 - 1024/9 PolyLog[2,1/2] Zeta[2])
   */

   const double cSM20 = 124.0607590216004;

   // Eq.(43) from arXiv:1807.03509
   const double dl3_analytic = alt2*as2*sb4*(
      1./27. * (
                6082 - 27832*LS + 14856*LS2 - 4032*LS3
                - 15408*z3 + 1728*LS*z3 - 27*cSM20
                + xt*(7616*LS - 11712*LS2 + 32*(-940 + 477*z3))
                + xt2*(28848 - 2640*LS + 1008*LS2 - 11880*z3)
                + xt3*(160*LS + 864*LS2 + 8*(2722 - 2259*z3))
                )
      );

   // Eq.(44) from arXiv:1807.03509
   const double dl3_analytic_shift = alt2*as2*sb4*(
      1./27. * (
                26916*LS - 18816*LS2 - 5904*LS3
                - xt*(-3744 + 14016*LS + 18816*LS2)
                - xt2*(29652 - 5424*LS - 9936*LS2)
                - xt3*(-6768 - 13152*LS - 2688*LS2)
                )
      );

   const int hierarchy = ho.getSuitableHierarchy();
   const auto hierarchy_name = ho.getH3mHierarchyNotation(hierarchy);

   INFO("MS = " << MS << ", tb = " << tb << ", xt = " << xt);
   INFO("yt = " << yt << ", g3 = " << g3 << ", LS = " << LS);
   INFO("αt^2 = " << alt2 << ", as^2 = " << as2);
   INFO("hierarchy = " << hierarchy_name);

   const double eps = 1e-5;

   CHECK_CLOSE(dl3, dl3_analytic, eps);
   CHECK_CLOSE(dl3_shift, dl3_analytic_shift, eps);
}

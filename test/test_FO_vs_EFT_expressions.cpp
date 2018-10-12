#include "doctest.h"
#include "Mh2EFTCalculator.hpp"
#include "MSSM_mass_eigenstates.hpp"
#include "EFTFlags.hpp"

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
   pars.Ad    << 0, 0, 0, 0 ,0 ,0 ,0 ,0, 0;
   pars.Yu    << 0, 0, 0, 0, 0, 0, 0, 0, 0.8;
   pars.Ye    << 0, 0, 0, 0, 0, 0, 0, 0, 0.1;
   pars.Yd    << 0, 0, 0, 0 ,0 ,0 ,0 ,0, 0.3;
   pars.M1    = MS;
   pars.M2    = MS;
   pars.MG    = mg;
   // pars.MW    = 74.597;
   // pars.MZ    = 85.7704;
   // pars.Mt    = 144.337;
   // pars.Mb    = 2.37054;
   // pars.Mtau  = 1.2;
   pars.MA    = MS;

   pars.validate(false);

   return pars;
}

himalaya::Parameters make_gaugeless(const himalaya::Parameters& pars)
{
   auto gl = pars;
   gl.g1 = 0.;
   gl.g2 = 0.;
   return gl;
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
   mhc.setCorrectionFlag(EFTOrders::G12G22, 0);
   mhc.setCorrectionFlag(EFTOrders::G12YB2, 0);
   mhc.setCorrectionFlag(EFTOrders::G14, 0);
   mhc.setCorrectionFlag(EFTOrders::G24, 0);
   mhc.setCorrectionFlag(EFTOrders::G12YB2, 0);
   mhc.setCorrectionFlag(EFTOrders::G22YB2, 0);
   mhc.setCorrectionFlag(EFTOrders::YB4, 0);
   mhc.setCorrectionFlag(EFTOrders::G12YTAU2, 0);
   mhc.setCorrectionFlag(EFTOrders::G22YTAU2, 0);
   mhc.setCorrectionFlag(EFTOrders::YTAU4, 0);
   mhc.setCorrectionFlag(EFTOrders::G12YT2, 0);
   mhc.setCorrectionFlag(EFTOrders::G22YT2, 0);

   return mhc.getDeltaMh2EFT1Loop(1,1);
}

} // anonymous namespace

TEST_CASE("test_FO_1loop_gaugeless")
{
   using namespace himalaya::mh1l;

   const double eps = 1e-10;
   const auto p    = make_point();
   const auto p_gl = make_gaugeless(p);
   const MSSM_mass_eigenstates me(p);
   const MSSM_mass_eigenstates me_gl(p_gl); // gaugeless
   const auto DMh2_1 = me.delta_mh2_1loop_gaugeless();
   const auto DMh2_2 = me_gl.delta_mh2_1loop_gaugeless();
   const auto DMh2_3 = me_gl.delta_mh2_1loop(0);

   CHECK_CLOSE(DMh2_1(0,0), DMh2_2(0,0), eps);
   CHECK_CLOSE(DMh2_1(0,1), DMh2_2(0,1), eps);
   CHECK_CLOSE(DMh2_1(1,0), DMh2_2(1,0), eps);
   CHECK_CLOSE(DMh2_1(1,1), DMh2_2(1,1), eps);
   CHECK_CLOSE(DMh2_2(0,0), DMh2_3(0,0), 1e-7);
   CHECK_CLOSE(DMh2_2(0,1), DMh2_3(0,1), 1e-7);
   CHECK_CLOSE(DMh2_2(1,0), DMh2_3(1,0), 1e-7);
   CHECK_CLOSE(DMh2_2(1,1), DMh2_3(1,1), 1e-7);
}

TEST_CASE("test_EFT_vs_FO_1loop")
{
   using namespace himalaya;
   using namespace himalaya::mh1l;
   using namespace himalaya::mh2_eft;

   const auto p = make_point();

   const auto Mh2_EFT_0L = calc_Mh2_EFT_0L(p);
   const auto Mh2_EFT_1L = calc_Mh2_EFT_1L(p);

   const MSSM_mass_eigenstates me(p);
   const auto Mh2_full_0L = me.calculate_Mh2(0);
   const auto Mh2_full_1L = me.calculate_Mh2(1);

   CHECK_CLOSE(Mh2_EFT_0L, Mh2_full_0L(0), 1e-6);
   // TODO increase precision
   CHECK_CLOSE(Mh2_EFT_1L, Mh2_full_1L(0), 1e-4);
}

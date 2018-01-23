#define BOOST_TEST_MODULE test_lambda_normalization

#include <boost/test/unit_test.hpp>
#include "HierarchyCalculator.hpp"
#include "Mh2EFTCalculator.hpp"

namespace {

double pow2(double x) { return x*x; }

himalaya::Parameters make_point()
{
   const double MS = 100000.;
   const double xt = 2.;
   const double tb = 5.;
   const double beta = std::atan(tb);
   const double v = 245.;

   const double eps = 0.001;

   himalaya::Parameters pars;
   pars.scale = MS;
   pars.mu    = MS;
   pars.g3    = 1.05733;
   pars.vu    = v*std::sin(beta);
   pars.vd    = v*std::cos(beta);
   pars.mq2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS)*(1 + eps);
   pars.md2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS);
   pars.mu2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS)*(1 - eps);
   pars.At    = xt*MS + pars.mu/tb;
   pars.Ab    = 0;
   pars.MG    = MS;
   pars.MW    = 74.597;
   pars.MZ    = 85.7704;
   pars.Mt    = 144.337;
   pars.Mb    = 2.37054;
   pars.MA    = MS;

   return pars;
}

// calculates a_t
double calc_at(const himalaya::Parameters& pars)
{
   const double Pi = 3.141592653589793;
   const double yt = std::sqrt(2.)*pars.Mt/pars.vu;
   const double tb = pars.vu/pars.vd;
   const double beta = std::atan(tb);

   return pow2(yt*std::sin(beta))/(4*Pi);
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
   const double at = calc_at(pars);
   const double tb = pars.vu/pars.vd;
   const double mt = pars.Mt;
   const double mQ32 = pars.mq2(2,2);
   const double mU32 = pars.mu2(2,2);
   const double Xt = pars.At - pars.mu/tb;
   const double MR2 = pow2(pars.scale);

   himalaya::mh2_eft::Mh2EFTCalculator mhc;

   return mhc.Mh2_EFT_1loop(at, mt, mQ32, mU32, Xt, MR2);
}

/// calculates Mh^2 in the EFT at 2-loop level
double calc_Mh2_EFT_2L(const himalaya::Parameters& pars)
{
   const double at = calc_at(pars);
   const double tb = pars.vu/pars.vd;
   const double mt = pars.Mt;
   const double mQ32 = pars.mq2(2,2);
   const double mU32 = pars.mu2(2,2);
   const double Xt = pars.At - pars.mu/tb;
   const double MR2 = pow2(pars.scale);
   const double g3 = pars.g3;
   const double m3 = pars.MG;

   himalaya::mh2_eft::Mh2EFTCalculator mhc;

   return mhc.Mh2_EFT_2loop(at, mt, mQ32, mU32, Xt, MR2, g3, m3);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_lambda_normalization)
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

   BOOST_TEST_MESSAGE("Mh(0L) = " << std::sqrt(Mh2_0L));
   BOOST_TEST_MESSAGE("Mh(1L) = " << std::sqrt(Mh2_1L));
   BOOST_TEST_MESSAGE("Mh(2L) = " << std::sqrt(Mh2_2L));
   BOOST_TEST_MESSAGE("Mh(3L) = " << std::sqrt(Mh2_3L));

   BOOST_CHECK_CLOSE_FRACTION(Mh2_0L, calc_Mh2_EFT_0L(pars), 1e-5);
   BOOST_CHECK_CLOSE_FRACTION(Mh2_1L, calc_Mh2_EFT_0L(pars)
                                    + calc_Mh2_EFT_1L(pars), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(Mh2_2L, calc_Mh2_EFT_0L(pars)
                                    + calc_Mh2_EFT_1L(pars)
                                    + calc_Mh2_EFT_2L(pars), 1e-6);

}

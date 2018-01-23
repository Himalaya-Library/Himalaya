#define BOOST_TEST_MODULE test_lambda_normalization

#include <boost/test/unit_test.hpp>
#include "HierarchyCalculator.hpp"

namespace {

double pow2(double x) { return x*x; }

himalaya::Parameters make_point()
{
   const double MS = 100000.;
   const double xt = 2.;
   const double tb = 5.;
   const double beta = std::atan(tb);
   const double v = 245.;

   himalaya::Parameters pars;
   pars.scale = MS;
   pars.mu    = MS;
   pars.g3    = 1.05733;
   pars.vu    = v*std::sin(beta);
   pars.vd    = v*std::cos(beta);
   pars.mq2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS);
   pars.md2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS);
   pars.mu2   << pow2(MS), 0, 0, 0, pow2(MS), 0, 0, 0, pow2(MS);
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

/// calculates Mh^2 in the decoupling limit
double calc_Mh2_tree_dec(const himalaya::Parameters& pars)
{
   const double beta = std::atan(pars.vu/pars.vd);
   const double cos_2beta = std::cos(2*beta);

   return pow2(pars.MZ * cos_2beta);
}

/// calculates lightest mass eigenvalue of given mass matrix
double calc_Mh2(const Eigen::Matrix2d& Mh_mat)
{
   const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(Mh_mat);
   return es.eigenvalues()(0);
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

   BOOST_CHECK_CLOSE_FRACTION(Mh2_0L, calc_Mh2_tree_dec(pars), 1e-5);
}

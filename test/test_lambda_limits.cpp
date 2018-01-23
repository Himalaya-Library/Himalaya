#define BOOST_TEST_MODULE test_lambda_limits

#include <boost/test/unit_test.hpp>
#include "HierarchyCalculator.hpp"

namespace {

himalaya::Parameters make_point(double eps = 0.)
{
   // MS = 1 TeV, Xt = -2 MS, TB = 5
   himalaya::Parameters pars;
   pars.scale = 1000;
   pars.mu    = 1000;
   pars.g3    = 1.05733;
   pars.vd    = 45.8309;
   pars.vu    = 229.155;
   pars.mq2   << 1e+06, 0, 0, 0, 1e+06, 0, 0, 0, 1e+06*(1 + eps);
   pars.md2   << 1e+06, 0, 0, 0, 1e+06, 0, 0, 0, 1e+06;
   pars.mu2   << 1e+06, 0, 0, 0, 1e+06, 0, 0, 0, 1e+06*(1 - eps);
   pars.At    = -1800;
   pars.Ab    = 5000;
   pars.MG    = 1000*(1 - eps);
   pars.MW    = 74.597;
   pars.MZ    = 85.7704;
   pars.Mt    = 144.337;
   pars.Mb    = 2.37054;
   pars.MA    = 1000;

   return pars;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_lambda_limit_degenerate)
{
   using namespace himalaya;

   const auto pars = make_point(0.03);
   auto hc = HierarchyCalculator(pars);
   const auto ho = hc.calculateDMh3L(false);

   const auto z2_deg          = ho.getZetaDegenerated();
   const auto z2_gen_Himalaya = ho.getZetaHimalaya();
   const auto z2_gen_EFT      = ho.getZetaEFT();

   BOOST_CHECK_CLOSE_FRACTION(z2_gen_Himalaya, z2_deg, 1e-5);
   BOOST_CHECK_CLOSE_FRACTION(z2_gen_EFT     , z2_deg, 1e-5);
}

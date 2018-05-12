#include "doctest.h"
#include "HierarchyCalculator.hpp"

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

namespace {

// point with small Xt, because of missing higher order Xt terms in H3m
himalaya::Parameters make_point(double eps = 0.)
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
   pars.At    = 300;
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

TEST_CASE("test_lambda_limit_degenerate")
{
   using namespace himalaya;

   const auto pars = make_point(0.03);
   auto hc = HierarchyCalculator(pars);
   const auto ho = hc.calculateDMh3L(false);

   const auto z2_gen_Himalaya = ho.getDeltaLambdaHimalaya();
   const auto z2_gen_EFT      = ho.getDeltaLambdaEFT();
   const auto uncertainty     = ho.getDeltaLambdaUncertainty();

   INFO("uncertainty: " << uncertainty);

   CHECK_CLOSE(z2_gen_Himalaya, z2_gen_EFT, uncertainty);
}

#include "doctest.h"
#include "HierarchyCalculator.hpp"
#include <cmath>

namespace {

himalaya::Parameters setup_SPS1a(){
   himalaya::Parameters pars;

   pars.scale = 4.67491329E+02;
   pars.mu = 3.52600579E+02;
   pars.g3 = 1.09949966E+00;
   pars.vd = 2.49832484E+01;
   pars.vu = 2.43650538E+02;
   pars.mq2 << 2.99793928E+05, 0, 0,
               0, 2.99792102E+05, 0,
               0, 0, 2.49327504E+05;
   pars.md2 << 2.78275669E+05, 0, 0,
               0, 2.78273780E+05, 0,
               0, 0, 2.74928741E+05;
   pars.mu2 << 2.80477426E+05, 0, 0,
               0, 2.80475621E+05, 0,
               0, 0, 1.80478484E+05;
   pars.Ad << 0, 0, 0, 0, 0, 0, 0, 0,  -798.8142296644842;
   pars.Au << 0, 0, 0, 0, 0, 0, 0, 0,  -506.4162662374052;

   pars.MA = 3.92960954E+02;
   pars.MG = 5.88220143E+02;
   pars.MW = 8.04136643E+01;
   pars.MZ = 9.06817306E+01;
   pars.Mt = 1.52117491E+02;
   pars.Mb = 2.42010269E+00;

   // let Himalaya calculate the stop masses and mixing:
   // pars.MSt << 3.83255677E+02, 5.70240743E+02;
   // pars.MSb << 4.98792475E+02, 5.28677648E+02;
   // pars.s2t = sin(2*asin(5.57720315E-01));
   // pars.s2b = sin(2*asin(-9.33860207E-01));

   return pars;
}

} // anonymous namespace

TEST_CASE("test_example")
{
   using RM22 = Eigen::Matrix2d;

   try {
      const himalaya::Parameters pars = setup_SPS1a();
      INFO(pars);

      himalaya::HierarchyCalculator hc(pars);
      const himalaya::HierarchyObject ho = hc.calculateDMh3L(false);

      const RM22 dMh0L = ho.getDMh(0);
      const RM22 dMh1L = ho.getDMh(1);
      const RM22 dMh2L = ho.getDMh(2);
      const RM22 dMh3L = ho.getDMh(3);
      const RM22 Mh2_mat = (dMh0L + dMh1L + dMh2L + dMh3L).eval();

      const Eigen::SelfAdjointEigenSolver<RM22> es(Mh2_mat);
      const auto Mh2 = es.eigenvalues();

      INFO("Mh^2 = " << Mh2);
      CHECK(std::isfinite(Mh2(0)));
      CHECK(std::isfinite(Mh2(1)));

      const double lam_3L_DRp = ho.getDLambdaEFT();
      const double Dlam_3L    = ho.getDLambdaEFTUncertainty();
      const double lam_3L_MS  = lam_3L_DRp + ho.getDLambdaEFTDRbarPrimeToMSbarShift();

      INFO("Δλ(3L,DR'-bar) = " << lam_3L_DRp << " +- " << Dlam_3L);
      INFO("Δλ(3L,MS-bar)  = " << lam_3L_MS  << " +- " << Dlam_3L);
      CHECK(std::isfinite(lam_3L_DRp));
      CHECK(std::isfinite(lam_3L_MS));
      CHECK(std::isfinite(Dlam_3L));
   } catch (const std::exception& e) {
      MESSAGE(e.what());
      REQUIRE_FALSE(false);
   }
}

himalaya::Parameters make_point(double MS, double xt, double tb)
{
   himalaya::Parameters pars;

   const double MS2 = MS*MS;
   const double Xt = xt*MS;
   const double beta = std::atan(tb);

   pars.scale = MS;
   pars.mu = MS;
   pars.g1 = 0.46;
   pars.g2 = 0.65;
   pars.g3 = 1.166;
   pars.vd = 246*std::cos(beta);
   pars.vu = 246*std::sin(beta);
   pars.mq2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.md2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.mu2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.ml2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.me2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.Au(2,2) = Xt + pars.mu/tb;
   pars.Yu(2,2) = 0.862;
   pars.Yd(2,2) = 0.133;
   pars.Ye(2,2) = 0.101;
   pars.MA = MS;
   pars.M1 = MS;
   pars.M2 = MS;
   pars.MG = MS;

   return pars;
}

TEST_CASE("test_Himalaya_2.0_example_1")
{
   // create parameter point
   const auto point = make_point(2000, std::sqrt(6.), 20);

   // create calculator object
   himalaya::HierarchyCalculator hc(point);

   // calculate 3-loop corrections
   const auto ho = hc.calculateDMh3L(false);

   // get 3-loop contribution to CP-even Higgs mass matrix
   const auto dMh_3L = ho.getDMh(3);

   // get 3-loop contribution to lambda
   const auto delta_lambda_3L = ho.getDLambdaEFT();

   // get uncertainty for 3-loop lambda
   const auto delta_lambda_3L_uncertainty = ho.getDLambdaEFTUncertainty();
}

#include "doctest.h"

#include "himalaya/mh2_eft/Mh2EFTCalculator.hpp"
#include "himalaya/mh2_eft/ThresholdCalculator.hpp"
#include "himalaya/misc/CouplingOrders.hpp"

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

namespace {
   const double Pi  = 3.1415926535897932384626433832795;
   template <typename T> T pow2(T x)  { return x*x; }
   template <typename T> T pow4(T x)  { return x*x*x*x; }
   template <typename T> T pow6(T x)  { return x*x*x*x*x*x; }

void _test_EFT_expressions(const himalaya::Parameters& p)
{
   using namespace himalaya;
   using namespace himalaya::mh2_eft;

   Mh2EFTCalculator mhc(p);
   ThresholdCalculator tc(p, true);
   const double eps = 1e-5;
   const double yt = sqrt(2)*p.Mt/std::sqrt(pow2(p.vu) + pow2(p.vd));

   {
      const double delta_mh_1l = mhc.getDeltaMh2EFT1Loop(1, 1);
      CHECK_CLOSE(delta_mh_1l, 4664.4945767867, eps);
   }

   {
      mhc.setCorrectionFlag(CouplingOrders::G32YB4, 0);
      mhc.setCorrectionFlag(CouplingOrders::YB6, 0);
      mhc.setCorrectionFlag(CouplingOrders::YT6, 0);
      mhc.setCorrectionFlag(CouplingOrders::YTAU2YB4, 0);
      mhc.setCorrectionFlag(CouplingOrders::YTAU4YB2, 0);
      mhc.setCorrectionFlag(CouplingOrders::YTAU6, 0);
      mhc.setCorrectionFlag(CouplingOrders::YT2YB4, 0);
      mhc.setCorrectionFlag(CouplingOrders::YB2YT4, 0);
      const double delta_mh_2l = mhc.getDeltaMh2EFT2Loop(1, 1);
      const double pref = 1./pow4(4*Pi) * pow2(p.Mt * yt * p.g3);
      CHECK_CLOSE(delta_mh_2l, pref*778.7287955, eps);
   }

   {
      const double delta_mh_3l = mhc.getDeltaMh2EFT3Loop(1, 1);
      const double susy_logs_3l = mhc.getDeltaMh2EFT3Loop(0, 1);
      CHECK_CLOSE(delta_mh_3l, 136.6811130551, eps);
      CHECK_CLOSE(susy_logs_3l, 11.7401841057, eps);
   }

   // check threshold corrections
   {

      const double g3as = tc.getThresholdCorrection(ThresholdCouplingOrders::G3_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(g3as, -0.1731100882, eps);

      // general mass case
      tc.setLimit(Limits::GENERAL);
      const double ytas_lim_gen = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_gen, -3.215927532, eps);
      const double ytas2_lim_gen = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_gen, -55.89911582, eps);
      const double lamat_lim_gen = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamat_lim_gen, 14.31721183, eps);
      const double lamatas_lim_gen = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_gen, -87.06409990, eps);
      const double lamatas2_lim_gen = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_gen, 207.7550117585, eps);
      const double drToMsShift_gen = tc.getDRbarPrimeToMSbarShift(3, 1);
      CHECK_CLOSE(drToMsShift_gen, -4573.1448163138, eps);

      // mQ3 = m3
      tc.setLimit(Limits::MQ3_EQ_M3);
      const double ytas_lim_mq3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_mq3_m3, -2.975633850, eps);
      const double ytas2_lim_mq3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_mq3_m3, -54.3868694883, eps);
      const double lamatas_lim_mq3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_mq3_m3, -71.71873062, eps);
      const double lamatas2_lim_mq3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_mq3_m3, 71.0964678625, eps);
      const double drToMsShift_mq3_m3 = tc.getDRbarPrimeToMSbarShift(3, 1);
      CHECK_CLOSE(drToMsShift_mq3_m3, -4475.2303977325, eps);

      // mQ3 = mU3
      tc.setLimit(Limits::MQ3_EQ_MU3);
      const double ytas_lim_mq3_mu3 = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_mq3_mu3, -3.187793781, eps);
      const double ytas2_lim_mq3_mu3 = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_mq3_mu3, -55.529811961, eps);
      const double lamat_lim_mq3_mu3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamat_lim_mq3_mu3, 14.35556635, eps);
      const double lamatas_lim_mq3_mu3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_mq3_mu3, -102.3684261, eps);
      const double lamatas2_lim_mq3_mu3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_mq3_mu3, 231.382248943, eps);
      const double drToMsShift_mq3_mu3 = tc.getDRbarPrimeToMSbarShift(3, 1);
      CHECK_CLOSE(drToMsShift_mq3_mu3, -4044.3655622086, eps);

      // mq3 = mu3 = m3
      tc.setLimit(Limits::MQ3_EQ_MU3_EQ_M3);
      const double ytas_lim_mq3_mu3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_mq3_mu3_m3, -2.957820456, eps);
      const double ytas2_lim_mq3_mu3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_mq3_mu3_m3, -53.982275575, eps);
      const double lamatas_lim_mq3_mu3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_mq3_mu3_m3, -87.55800303, eps);
      const double lamatas2_lim_mq3_mu3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_mq3_mu3_m3, 102.8713758798, eps);
      const double drToMsShift_mq3_mu3_m3 = tc.getDRbarPrimeToMSbarShift(3, 1);
      CHECK_CLOSE(drToMsShift_mq3_mu3_m3, -3974.9639012363, eps);

      // mU3 = m3
      tc.setLimit(Limits::MU3_EQ_M3);
      const double ytas_lim_mu3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_mu3_m3, -2.747511273, eps);
      const double ytas2_lim_mu3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_mu3_m3, -52.7089852536, eps);
      const double lamatas_lim_mu3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_mu3_m3, -61.18314361, eps);
      const double lamatas2_lim_mu3_m3 = tc.getThresholdCorrection(ThresholdCouplingOrders::LAMBDA_AT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_mu3_m3, -32.4017074174, eps);
      const double drToMsShift_mu3_m3 = tc.getDRbarPrimeToMSbarShift(3, 1);
      CHECK_CLOSE(drToMsShift_mu3_m3, -4323.0510347014, eps);

      // degenerate mass case
      tc.setLimit(Limits::DEGENERATE);
      const double ytas2_lim_degen = tc.getThresholdCorrection(ThresholdCouplingOrders::YT_AS2,
	 RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_degen, -52.81677866, eps);
      const double drToMsShift_degen = tc.getDRbarPrimeToMSbarShift(3, 1);
      CHECK_CLOSE(drToMsShift_degen, -3900.7300573389, eps);
   }
}

} // anonymous namespace

himalaya::Parameters test_point(){
   himalaya::Parameters pars;

   pars.scale = 4.67491329E+02;
   pars.mu = 3.52600579E+02;
   pars.g1 = 0.441854;
   pars.g2 = 0.656631;
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
   pars.ml2 << 2.99793928E+05, 0, 0,
               0, 2.99792102E+05, 0,
               0, 0, 2.49327504E+05;
   pars.me2 << 2.99793928E+05, 0, 0,
               0, 2.99792102E+05, 0,
               0, 0, 2.49327504E+05;

   pars.Ad << 0, 0, 0, 0, 0, 0, 0, 0,  -798.8142296644842;
   pars.Au << 0, 0, 0, 0, 0, 0, 0, 0,  -506.4162662374052;
   pars.Ae << 0, 0, 0, 0, 0, 0, 0, 0,   -850.;

   pars.M1 = 500.;
   pars.M2 = 600.;
   pars.MA = 3.92960954E+02;
   pars.MG = 5.88220143E+02;
   pars.MW = 8.04136643E+01;
   pars.MZ = 9.06817306E+01;
   pars.Mt = 1.52117491E+02;
   pars.Mb = 2.42010269E+00;
   pars.Mtau = 1.2;
   pars.MSt << 3.83255677E+02, 5.70240743E+02;
   pars.MSb << 4.98792475E+02, 5.28677648E+02;
   pars.s2t = sin(2*asin(5.57720315E-01));
   pars.s2b = sin(2*asin(-9.33860207E-01));

   pars.massLimit3LThreshold = static_cast<int>(himalaya::mh2_eft::Limits::GENERAL);

   return pars;
}

TEST_CASE("test_EFT_expressions")
{
   const himalaya::Parameters p = test_point();
    _test_EFT_expressions(p);
}

#include "doctest.h"

#define private public
#include "Mh2EFTCalculator.hpp"
#include "ThresholdCalculator.hpp"
#include "Hierarchies.hpp"
#include <iostream>
#include <iomanip>
#undef private

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

namespace {

void _test_EFT_expressions(const himalaya::Parameters& p, double msq2)
{
   himalaya::mh2_eft::Mh2EFTCalculator mhc(p, msq2);
   himalaya::ThresholdCalculator tc(p, msq2, true);
   const double eps = 1e-5;
   
   {
      const double delta_mh_1l = mhc.getDeltaMh2EFT1Loop(1, 1);
      const double susy_logs_1l = mhc.getDeltaMh2EFT1Loop(0, 1);
      CHECK_CLOSE(delta_mh_1l, 41.26267451, eps);
      CHECK_CLOSE(susy_logs_1l, 14.31721183, eps);
   }
   
   {
      const double delta_mh_2l = mhc.getDeltaMh2EFT2Loop(1, 1);
      const double susy_logs_2l = mhc.getDeltaMh2EFT2Loop(0, 1);
      CHECK_CLOSE(delta_mh_2l, 778.7287955, eps);
      CHECK_CLOSE(susy_logs_2l, 19.92610215, eps);
   }
   
   {
      const double delta_mh_3l = mhc.getDeltaMh2EFT3Loop(1, 1);
      const double susy_logs_3l = mhc.getDeltaMh2EFT3Loop(0, 1);
      CHECK_CLOSE(delta_mh_3l, 20632.69911, eps);
      CHECK_CLOSE(susy_logs_3l, 1772.341613, eps);
   }
   
   // check threshold corrections
   {
      const double Xt4 = 8.666148873*1e10;
      const double Xt5 = -4.702000616*1e13;
      const double Xt6 = 2.551168935*1e16;
      
      const double g3as = tc.getThresholdCorrection(himalaya::ThresholdVariables::G3_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(g3as, -0.1736769306, eps);
      
      // general mass case
      tc.setLimit(himalaya::Limits::GENERAL);
      const double ytas_lim_gen = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_gen, -3.215927532, eps);
      const double ytas2_lim_gen = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_gen, -55.90563134, eps);
      const double lamat_lim_gen = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamat_lim_gen, 14.31721183, eps);
      const double lamatas_lim_gen = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_gen, -87.06409990, eps);
      const double lamatas2_lim_gen = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_gen, 207.7740287, eps);
      const double xt4_lim_gen = tc.getXtTerms(4, 1)*1e8/Xt4;
      CHECK_CLOSE(xt4_lim_gen, -1.035236728, eps);
      const double xt5_lim_gen = tc.getXtTerms(5, 1)*1e12/Xt5;
      CHECK_CLOSE(xt5_lim_gen, 2.741433435, eps);
      const double xt6_lim_gen = tc.getXtTerms(6, 1)*1e15/Xt6;
      CHECK_CLOSE(xt6_lim_gen, 2.094602505, eps);
      
      // mQ3 = m3
      tc.setLimit(himalaya::Limits::MQ3_EQ_M3);
      const double ytas_lim_mq3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_mq3_m3, -2.975633850, eps);
      const double ytas2_lim_mq3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_mq3_m3, -54.39344951, eps);
      const double lamatas_lim_mq3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_mq3_m3, -71.71873062, eps);
      const double lamatas2_lim_mq3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_mq3_m3, 71.11240105, eps);
      const double xt4_lim_mq3_m3 = tc.getXtTerms(4, 1)*1e9/Xt4;
      CHECK_CLOSE(xt4_lim_mq3_m3, -9.723929768, eps);
      const double xt5_lim_mq3_m3 = tc.getXtTerms(5, 1)*1e12/Xt5;
      CHECK_CLOSE(xt5_lim_mq3_m3, 2.552971369, eps);
      const double xt6_lim_mq3_m3 = tc.getXtTerms(6, 1)*1e15/Xt6;
      CHECK_CLOSE(xt6_lim_mq3_m3, 1.922855926, eps);
      
      // mQ3 = mU3
      tc.setLimit(himalaya::Limits::MQ3_EQ_MU3);
      const double ytas_lim_mq3_mu3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_mq3_mu3, -3.187793781, eps);
      const double ytas2_lim_mq3_mu3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_mq3_mu3, -55.53609836, eps);
      const double lamat_lim_mq3_mu3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamat_lim_mq3_mu3, 14.35556635, eps);
      const double lamatas_lim_mq3_mu3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_mq3_mu3, -102.3684261, eps);
      const double lamatas2_lim_mq3_mu3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_mq3_mu3, 231.4042763, eps);
      const double xt4_lim_mq3_mu3 = tc.getXtTerms(4, 1)*1e9/Xt4;
      CHECK_CLOSE(xt4_lim_mq3_mu3, -7.294674107, eps);
      const double xt5_lim_mq3_mu3 = tc.getXtTerms(5, 1)*1e12/Xt5;
      CHECK_CLOSE(xt5_lim_mq3_mu3, 1.595070724, eps);
      const double xt6_lim_mq3_mu3 = tc.getXtTerms(6, 1)*1e15/Xt6;
      CHECK_CLOSE(xt6_lim_mq3_mu3, 1.264126661, eps);
      
      // mq3 = mu3 = m3
      tc.setLimit(himalaya::Limits::MQ3_EQ_MU3_EQ_M3);
      const double ytas_lim_mq3_mu3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_mq3_mu3_m3, -2.957820456, eps);
      const double ytas2_lim_mq3_mu3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_mq3_mu3_m3, -53.98859513, eps);
      const double lamatas_lim_mq3_mu3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_mq3_mu3_m3, -87.55800303, eps);
      const double lamatas2_lim_mq3_mu3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_mq3_mu3_m3, 102.8902455, eps);
      const double xt4_lim_mq3_mu3_m3 = tc.getXtTerms(4, 1)*1e9/Xt4;
      CHECK_CLOSE(xt4_lim_mq3_mu3_m3, -6.800716374, eps);
      const double xt5_lim_mq3_mu3_m3 = tc.getXtTerms(5, 1)*1e12/Xt5;
      CHECK_CLOSE(xt5_lim_mq3_mu3_m3, 1.509254564, eps);
      const double xt6_lim_mq3_mu3_m3 = tc.getXtTerms(6, 1)*1e15/Xt6;
      CHECK_CLOSE(xt6_lim_mq3_mu3_m3, 1.147009210, eps);
      
      // mU3 = m3
      tc.setLimit(himalaya::Limits::MU3_EQ_M3);
      const double ytas_lim_mu3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas_lim_mu3_m3, -2.747511273, eps);
      const double ytas2_lim_mu3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_mu3_m3, -52.71555755, eps);
      const double lamatas_lim_mu3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas_lim_mu3_m3, -61.18314361, eps);
      const double lamatas2_lim_mu3_m3 = tc.getThresholdCorrection(himalaya::ThresholdVariables::LAMBDA_AT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(lamatas2_lim_mu3_m3, -32.38839351, eps);
      const double xt4_lim_mu3_m3 = tc.getXtTerms(4, 1)*1e9/Xt4;
      CHECK_CLOSE(xt4_lim_mu3_m3, -8.949769828, eps);
      const double xt5_lim_mu3_m3 = tc.getXtTerms(5, 1)*1e12/Xt5;
      CHECK_CLOSE(xt5_lim_mu3_m3, 2.854004720, eps);
      const double xt6_lim_mu3_m3 = tc.getXtTerms(6, 1)*1e15/Xt6;
      CHECK_CLOSE(xt6_lim_mu3_m3, 1.726207002, eps);
      
      // degenerate mass case
      tc.setLimit(himalaya::Limits::DEGENERATE);
      const double ytas2_lim_degen = tc.getThresholdCorrection(himalaya::ThresholdVariables::YT_AS2,
	 himalaya::RenSchemes::TEST, 1);
      CHECK_CLOSE(ytas2_lim_degen, -52.81677866, eps);
      const double xt4_lim_degen = tc.getXtTerms(4, 1)*1e9/Xt4;
      CHECK_CLOSE(xt4_lim_degen, -6.758025736, eps);
      const double xt5_lim_degen = tc.getXtTerms(5, 1)*1e12/Xt5;
      CHECK_CLOSE(xt5_lim_degen, 1.470013998, eps);
   }
}

} // anonymous namespace

himalaya::Parameters test_point(){
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
   pars.Ab = -798.8142296644842;
   pars.At = -506.4162662374052;

   pars.MA = 3.92960954E+02;
   pars.MG = 5.88220143E+02;
   pars.MW = 8.04136643E+01;
   pars.MZ = 9.06817306E+01;
   pars.Mt = 1.52117491E+02;
   pars.Mb = 2.42010269E+00;
   pars.MSt << 3.83255677E+02, 5.70240743E+02;
   pars.MSb << 4.98792475E+02, 5.28677648E+02;
   pars.s2t = sin(2*asin(5.57720315E-01));
   pars.s2b = sin(2*asin(-9.33860207E-01));

   pars.massLimit3LThreshold = himalaya::Limits::GENERAL;
   
   return pars;
}

TEST_CASE("test_EFT_expressions")
{
   const himalaya::Parameters p = test_point();
   
   const double msq = 533.204;
   
    _test_EFT_expressions(p, msq * msq);
}

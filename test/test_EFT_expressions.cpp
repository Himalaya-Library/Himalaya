#define BOOST_TEST_MODULE test_EFT_expressions

#include <boost/test/unit_test.hpp>
#include "HierarchyCalculator.hpp"

#define private public
#include "Mh2EFTCalculator.hpp"
#undef private

namespace {

void _test_EFT_logs(double mQ32, double mU32, double Xt, double MR2, double m3, double msq2)
{
   himalaya::mh2_eft::Mh2EFTCalculator mhc;

   {
      const double c_as_0_log_0 = mhc.coeff_as_0_log_0(mQ32, mU32, Xt, MR2);
      BOOST_CHECK_CLOSE_FRACTION(c_as_0_log_0, 65.75234703, 1e-3);
   }

   {
      const double c_as_0_log_1 = mhc.coeff_as_0_log_1();
      BOOST_CHECK_CLOSE_FRACTION(c_as_0_log_1, 12., 1e-3);
   }

   {
      const double c_as_1_log_0 = mhc.coeff_as_1_log_0(mQ32, mU32, Xt, m3, MR2);
      BOOST_CHECK_CLOSE_FRACTION(c_as_1_log_0, 977.1683309, 1e-3);
   }

   {
      const double c_as_1_log_1 = mhc.coeff_as_1_log_1(mQ32, mU32, Xt, m3, MR2);
      BOOST_CHECK_CLOSE_FRACTION(c_as_1_log_1, 186.2109432, 1e-3);
   }

   {
      const double c_as_1_log_2 = mhc.coeff_as_1_log_2();
      BOOST_CHECK_CLOSE_FRACTION(c_as_1_log_2, 96., 1e-3);
   }

   {
      const double c_as_2_log_0 = mhc.coeff_as_2_log_0(mQ32, mU32, Xt, m3, msq2, MR2);
      BOOST_CHECK_CLOSE_FRACTION(c_as_2_log_0, -14454.64603, 1e-3);
   }

   {
      const double c_as_2_log_1 = mhc.coeff_as_2_log_1(mQ32, mU32, Xt, m3, msq2, MR2);
      BOOST_CHECK_CLOSE_FRACTION(c_as_2_log_1, -5680.883822, 1e-3);
   }

   {
      const double c_as_2_log_2 = mhc.coeff_as_2_log_2(mQ32, mU32, Xt, m3, msq2, MR2);
      BOOST_CHECK_CLOSE_FRACTION(c_as_2_log_2, 2877.986080, 1e-3);
   }

   {
      const double c_as_2_log_3 = mhc.coeff_as_2_log_3();
      BOOST_CHECK_CLOSE_FRACTION(c_as_2_log_3, 736., 1e-3);
   }
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_EFT_logs)
{
   const double mQ32 = 10000;
   const double mU32 = 20000;
   const double Xt = 2 * 100.;
   const double MR2 = 500;
   const double m3 = 300;
   const double msq2 = 400;

   _test_EFT_logs(mQ32, mU32, Xt, MR2, m3, msq2);
}

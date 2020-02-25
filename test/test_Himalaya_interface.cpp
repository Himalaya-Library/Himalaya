#include "doctest.h"

#include "himalaya/Himalaya_interface.hpp"

#define CHECK_CLOSE(a,b,eps) CHECK((a) == doctest::Approx(b).epsilon(eps))

namespace {

himalaya::Parameters test_stop_mixing() {
   himalaya::Parameters pars;

   const double MS = 10000.;
   const double MS2 = MS*MS;
   const double Xt = -sqrt(6.)*MS;
   const double tb = 20.;

   pars.scale = MS;
   pars.mu = MS;
   pars.g3 = 1.10073;
   pars.vd = 246*std::cos(std::atan(tb));
   pars.vu = 246*std::sin(std::atan(tb));
   pars.mq2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.md2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.mu2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, 9.99995e+07;
   pars.Ad(2,2) = 0;
   pars.Au(2,2) = Xt + pars.mu/tb;
   pars.MA = MS;
   pars.MG = 4*MS;
   pars.MW = 78.9441;
   pars.MZ = 90.5119;
   pars.Mt = 154.682;
   pars.Mb = 2.50901;

   return pars;
}

} // anonymous namespace

// This test ensures that for test_stop_mixing() point the sign of s2t
// is positive.  This must be the case, because mq2(2,2) < mu(2,2), so
// ~t1 = ~tL and ~t2 = ~tR.
TEST_CASE("test_stop_mixing")
{
   auto point = test_stop_mixing();
   point.validate(false);

   CHECK_CLOSE(point.s2t, 1., 0.001);
   CHECK(point.MSt(0) < point.MSt(1));
}

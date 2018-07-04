// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "HierarchyCalculator.hpp"
#include "Logger.hpp"

himalaya::Parameters setup_point(double MS, double tb, double xt)
{
   himalaya::Parameters pars;

   const double MS2 = MS*MS;
   const double Xt = xt*MS;
   const double beta = std::atan(tb);

   pars.scale = MS;
   pars.mu = MS;
   pars.g1 = 0.46;
   pars.g2 = 0.65;
   pars.g3 = 1.10073;
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
   pars.Au << 0, 0, 0,
              0, 0, 0, 0,
              0, Xt + pars.mu/tb;
   pars.Ad << 0, 0, 0,
              0, 0, 0,
              0, 0, 0;
   pars.Ae << 0, 0, 0,
              0, 0, 0,
              0, 0, 0;
   pars.MA = MS;
   pars.M1 = MS;
   pars.M2 = MS;
   pars.MG = MS;
   pars.Mt = 154.682;
   pars.Mb = 2.50901;
   pars.Mtau = 1.777;

   return pars;
}

int main()
{
   const std::vector<himalaya::Parameters> points = {
      setup_point(2000., 10., 0.1)
   };

   for (const auto& point: points) {
      try {
	 // init hierarchy calculator
	 himalaya::HierarchyCalculator hc(point);

	 // calculate the 3-loop corrections O(α_t*α_s^2)
	 const auto hoTop = hc.calculateDMh3L(false);

	 std::cout << hoTop;

	 // calculate the 3-loop corrections O(α_b*α_s^2)
	 // himalaya::HierarchyObject hoBot = hc.calculateDMh3L(true);
      } catch (const std::exception& e) {
         ERROR_MSG(e.what());
      }
   }

   return 0;
}

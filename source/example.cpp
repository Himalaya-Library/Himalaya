// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "HierarchyCalculator.hpp"
#include "Mh2EFTCalculator.hpp"
#include "Logger.hpp"
#include  <iomanip>

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
   pars.mq2 << MS2*1.1, 0, 0,
               0, MS2*1.2, 0,
               0, 0, MS2*1.3;
   pars.md2 << MS2*0.9, 0, 0,
               0, MS2*0.8, 0,
               0, 0, MS2*0.7;
   pars.mu2 << MS2*2, 0, 0,
               0, MS2*2.1, 0,
               0, 0, MS2*2.2;
   pars.ml2 << MS2*3, 0, 0,
               0, MS2*3.1, 0,
               0, 0, MS2*3.2;
   pars.me2 << MS2*4, 0, 0,
               0, MS2*4.1, 0,
               0, 0, MS2*4.2;
   pars.Au << 0, 0, 0,
              0, 0, 0, 0,
              0, Xt + pars.mu/tb;
   pars.Ad << 0, 0, 0,
              0, 0, 0,
              0, 0, 0;
   pars.Ae << 0, 0, 0,
              0, 0, 0,
              0, 0, 0;
   pars.MA = MS*3;
   pars.M1 = MS*5;
   pars.M2 = MS*6;
   pars.MG = MS*7;
   pars.Mt = 154.682;
   pars.Mb = 2.50901;
   pars.Mtau = 1.777;

   return pars;
}

int main()
{
   const std::vector<himalaya::Parameters> points = {
      setup_point(2500., 10., 0.)
   };

   for (const auto& point: points) {
      try {
         // init hierarchy calculator
         himalaya::HierarchyCalculator hc(point, 0);

         // calculate the 3-loop corrections O(α_t*α_s^2)
         //const auto hoTop = hc.calculateDMh3L(false);

         //std::cout << hoTop;
         // calculate the 3-loop corrections O(α_b*α_s^2)
         //himalaya::HierarchyObject hoBot = hc.calculateDMh3L(true);

         // check EFT expressions
         /*std::cout << std::setprecision(16) << point;
         himalaya::mh2_eft::Mh2EFTCalculator mh2EFTCalculator(point);
         std::cout << mh2EFTCalculator.getDeltaMh2EFT2Loop(1,1) << "\n";*/
      } catch (const std::exception& e) {
         ERROR_MSG(e.what());
      }
   }

   return 0;
}

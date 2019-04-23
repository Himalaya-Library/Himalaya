// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#include "HierarchyCalculator.hpp"
#include "Mh2EFTCalculator.hpp"
#include "MSSM_mass_eigenstates.hpp"
#include  <iostream>

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
   pars.Au << 0, 0, 0,
              0, 0, 0,
              0, 0, Xt + pars.mu/tb;
   pars.Ad << 0, 0, 0,
              0, 0, 0,
              0, 0, 0;
   pars.Ae << 0, 0, 0,
              0, 0, 0,
              0, 0, 0;
   pars.Yu << 0, 0, 0, 0, 0, 0, 0, 0, 0.862;
   pars.Yd << 0, 0, 0, 0 ,0 ,0 ,0 ,0, 0.133;
   pars.Ye << 0, 0, 0, 0, 0, 0, 0, 0, 0.101;
   pars.MA = MS;
   pars.M1 = MS;
   pars.M2 = MS;
   pars.MG = MS;

   pars.validate(true);

   return pars;
}

himalaya::Parameters make_gaugeless(const himalaya::Parameters& pars)
{
   auto gl = pars;
   gl.g1 = 0.;
   gl.g2 = 0.;
   gl.MW = himalaya::NaN;
   gl.MZ = himalaya::NaN;
   gl.validate(false);
   return gl;
}

int main()
{
   const std::vector<himalaya::Parameters> points = {
      setup_point(2000., 20., std::sqrt(6.))
   };

   for (const auto& point: points) {
      try {
         // init hierarchy calculator
         himalaya::HierarchyCalculator hc(point, 0);

         // calculate the 3-loop corrections O(α_t*α_s^2)
         const auto hoTop = hc.calculateDMh3L(false);
         std::cout << hoTop;

         // calculate the 3-loop corrections O(α_b*α_s^2)
         //himalaya::HierarchyObject hoBot = hc.calculateDMh3L(true);

         /*const auto point_gl = make_gaugeless(point);

         // calculate fixed-order corrections for v^2 << MS^2
         himalaya::mh2_eft::Mh2EFTCalculator meft(point_gl);
         //himalaya::mh2_eft::Mh2EFTCalculator meft(point);
         const auto dmh2_eft_0l = meft.getDeltaMh2EFT0Loop();
         const auto dmh2_eft_1l = meft.getDeltaMh2EFT1Loop(1,1);
         const auto dmh2_eft_2l = meft.getDeltaMh2EFT2Loop(1,1);

         std::cout << "Mh^2_EFT_0L  = " << dmh2_eft_0l << " GeV^2 O(g1^2, g2^2)\n";
         std::cout << "ΔMh^2_EFT_1L = " << dmh2_eft_1l << " GeV^2 O(αt + αb + ατ)\n";
         std::cout << "ΔMh^2_EFT_2L = " << dmh2_eft_2l
                   << " GeV^2 O((αt+ab)*αs + (αt+αb)^2 + ab*aτ + aτ^2)\n";

         // calculate fixed-order corrections
         himalaya::mh2_fo::MSSM_mass_eigenstates mfo(point_gl, true);
         //himalaya::mh2_fo::MSSM_mass_eigenstates mfo(point);
         const auto dmh_fo     = mfo.calculate_Mh2();
         const auto dmh2_fo_0l = std::get<0>(dmh_fo);
         const auto dmh2_fo_1l = std::get<1>(dmh_fo);
         const auto dmh2_fo_2l = std::get<2>(dmh_fo);

         std::cout << "Mh^2_FO_0L   = " << dmh2_fo_0l << " GeV^2 O(full)\n";
         std::cout << "ΔMh^2_FO_1L  = " << dmh2_fo_1l << " GeV^2 O(αt + αb + ατ)\n";
         std::cout << "ΔMh^2_FO_2L  = " << dmh2_fo_2l
                   << " GeV^2 O((αt+ab)*αs + (αt+αb)^2 + ab*aτ + aτ^2)\n";*/
      } catch (const std::exception& e) {
         std::cerr << "Error: " << e.what() << '\n';
      }
   }

   return 0;
}

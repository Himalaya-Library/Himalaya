#include "HierarchyCalculator.hpp"
#include "Himalaya_interface.hpp"
#include "HierarchyObject.hpp"
#include "Hierarchies.hpp"
#include "Logger.hpp"
#include <iostream>

himalaya::Parameters setup_point() {
   himalaya::Parameters pars;

   const double MS = 2000.;
   const double MS2 = MS*MS;
   const double Xt = 0.1*MS;
   const double tb = 10.;

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
               0, 0, MS2*0.9;
   pars.Ad << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   pars.Au << 0, 0, 0, 0, 0, 0, 0, 0, Xt + pars.mu/tb;
   pars.MA = MS;
   pars.MG = MS;
   pars.MW = 78.9441;
   pars.MZ = 90.5119;
   pars.Mt = 154.682;
   pars.Mb = 2.50901;

   return pars;
}

int main() {
   try{
      const std::vector<himalaya::Parameters> points = {
	 setup_point()
      }; 
      for (const auto& point: points) {
	 // init hierarchy calculator
	 himalaya::HierarchyCalculator hierarchyCalculator(point);

	 // calculate the 3-loop corrections with the suiatble hierarchy
	 // top
	 himalaya::HierarchyObject hoTop = hierarchyCalculator.calculateDMh3L(false);

	 std::cout << hoTop << "\n";
   
	 // bottom
	 //himalaya::HierarchyObject hoBot = hierarchyCalculator.calculateDMh3L(true);
      }
   }
   catch (const std::exception& e){
      ERROR_MSG(e.what());
      return EXIT_FAILURE;
   }
   return 0;
}

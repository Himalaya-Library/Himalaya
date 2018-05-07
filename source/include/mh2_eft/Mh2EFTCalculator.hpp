#pragma once

#include "Himalaya_interface.hpp"
#include "dilog.h"
#include <cmath>
//#include <limits>

namespace himalaya{
namespace mh2_eft{
   
   /**
    * The Mh2 EFT calculator class
    */
   class Mh2EFTCalculator{
   public:
      /**
       *	Constructor
       * 	@param p a HimalayaInterface struct
       * 	@param msq2 the averaged squark mass of the first two generations squared
       */
      Mh2EFTCalculator(const Parameters& p, const double msq2);
      /**
       * 	Constructor
       */
      Mh2EFTCalculator();
      /**
       * 	Returns the 1-loop EFT contribution to the Higgs mass
       * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
       * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
       */
      double getDeltaMh2EFT1Loop(int omitSMLogs, int omitMSSMLogs);
      /**
       * 	Returns the 1-loop EFT contribution to the Higgs mass
       * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
       * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
       */
      double getDeltaMh2EFT2Loop(int omitSMLogs, int omitMSSMLogs);
      /**
       * 	Returns the 1-loop EFT contribution to the Higgs mass
       * 	@param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
       * 	@param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
       * 	@param xtOrder an integer to subtract all Xt contributions starting at Xt^4 up to Xt^6 (e.g. xtOrder = 4 will subtract O(Xt^5 + Xt^6) terms). By default no Xt terms are subtracted. SM logs are set to 0 by default.
       */
      double getDeltaMh2EFT3Loop(int omitSMLogs, int omitMSSMLogs, int xtOrder = 7);
      /**
       *   Returns the matching relation of zeta_lambda^(2) for the degenerate mass case
       *   @param scale the renormalization scale
       *   @param mst1 the mass of the light stop quark
       *   @return zeta_lambda^(2)
       */
     double getZetaDegenerate(double scale, double mst1, double Xt, int omitlogs) const;
   private:
      /**
       * 	Returns a string of "true" or "false" corresponding to the argument
       * 	@param tf a boolean
       */
      std::string tf(const bool tf);
      Parameters p{};
      double msq2{};
   };
}	// mh2_eft
}	// himalaya

#pragma once

#include "Himalaya_interface.hpp"
#include "dilog.h"
#include <cmath>
#include <limits>
#include <map>

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
       * 	@param verbose a bool enable the output of the parameter validation. Enabled by default
       */
      Mh2EFTCalculator(const Parameters& p, double msq2 = std::numeric_limits<double>::quiet_NaN(), bool verbose = true);
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
       * 	@param omitDeltaLambda3L an integer flag to disable the MSSM contribution to delta_lambda_3L
       */
      double getDeltaMh2EFT3Loop(int omitSMLogs, int omitMSSMLogs, int omitDeltaLambda3L = 1);
      /**
       *   Returns the matching relation of delta_lambda 3L for the degenerate mass case
       *   @param scale the renormalization scale
       *   @param mst1 the mass of the light stop quark
       *   @return delta_Lambda 3L
       */
      double getDeltaLambdaDegenerate(double scale, double mst1, double Xt, int omitlogs) const;
      /**
       * 	Sets the flag to enable or disable a correction of a given variable
       * 	@param variable an integer taken from the EFTOrders enum
       * 	@param enable set to 1 to enable and to 0 to disable the chosen correction
       */
      void setCorrectionFlag(int variable, int enable);
   private:
      /**
       * 	Returns a string of "true" or "false" corresponding to the argument
       * 	@param tf a boolean
       */
      std::string tf(const bool tf);
      Parameters p{};
      double msq2{std::numeric_limits<double>::quiet_NaN()};
      std::map<unsigned int, unsigned int> orderMap{};	/** A map which holds all EFTOrders key value pairs to enable/disable certain corrections. */
   };
}	// mh2_eft
}	// himalaya

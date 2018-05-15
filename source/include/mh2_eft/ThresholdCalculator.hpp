#pragma once

#include "Himalaya_interface.hpp"
#include <limits>

namespace himalaya{
   
   class ThresholdCalculator{
   public:
      /**
       * 	Constructor
       * 	@param p a HimalayaInterface struct
       * 	@param msq2 the averaged squark mass of the first two generations squared
       * 	@param verbose a bool enable the output of the parameter validation. Enabled by default
       * 	@param check a boolean which indicates if the threshold corrections should be tested
       */
      ThresholdCalculator(const Parameters& p, double msq2 = std::numeric_limits<double>::quiet_NaN(), bool verbose = true, bool check = false);
      /**
       * 	Returns a specific threshold corrections for a given mass limit
       * 	@param variable an integer key for a threshold correctionn
       * 	@param scheme an integer key to set the scheme. Choices are {MSbar, DRbar', DRbar}
       * 	@param omitLogs an integer to omit all log mu terms
       * 	@return a threshold correction for a given variable in a given scheme for a suitable mass limit
       */
      double getThresholdCorrection(int variable, int scheme, int omitLogs);
      /**
       * 	Returns the O(Xt^n) contribution to Mh^2 without any prefactors at 3-loop level
       * 	@param xtOrder an integer key to get the Xt contributions at the given order starting at Xt^4 for xtOrder = 4
       * 	@param omitLogs an integer key to omit all log mu terms
       */
      double getXtTerms(int xtOrder, int omitLogs);
      /**
       * 	Returns the shift needed to convert the 3L threshold correction of lambda to the MSbar scheme
       * 	@param xtOrder an integer key to omit the Xt contributions starting at xtOrder + 1
       * 	@param omitLogs an integer key to omit all log mu terms
       * 	@param omitXtLogs an integer key to omit all Xt^4*Log[mu] and Xt^5*Log[mu] terms
       */
      double getDRbarPrimeToMSbarShift(int xtOrder, int omitLogs, int omitXtLogs = 1);
      /**
       * 	Sets the mass limit to check terms
       * 	@param limit an integer key for a mass limit
       */
      void setLimit(int limit);
   private:
      /**
       * 	Returns delta g3_as in the MSbar scheme for a given mass limit
       * 	@param omitLogs an integer key to omit all mu terms
       * 	@return delta g3_as in the MSbar scheme for a given mass limit
       */
      double getDeltaG3Alphas(int omitLogs);
      /**
       * 	Returns delta yt_as in the MSbar scheme for a given mass limit
       * 	@param limit an integer key for a mass limit
       * 	@param omitLogs an integer key to omit all mu terms
       * 	@return delta yt_as in the MSbar scheme for a given mass limit
       */
      double getDeltaYtAlphas(int limit, int omitLogs);
      /**
       * 	Returns delta yt_as^2 in the MSbar scheme for a given mass limit
       * 	@param limit an integer key for a mass limit
       * 	@param omitLogs an integer key to omit all mu terms
       * 	@return delta yt_as^2 in the MSbar scheme for a given mass limit
       */
      double getDeltaYtAlphas2(int limit, int omitLogs);
      /**
       * 	Returns delta lambda_at in the MSbar scheme for a given mass limit
       * 	@param limit an integer key for a mass limit
       * 	@param omitLogs an integer key to omit all mu terms
       * 	@return delta lambda_at in the MSbar scheme for a given mass limit
       */
      double getDeltaLambdaAlphat(int limit, int omitLogs);
      /**
       * 	Returns delta lambda_atas in the MSbar scheme for a given mass limit
       * 	@param limit an integer key for a mass limit
       * 	@param omitLogs an integer key to omit all mu terms
       * 	@return delta lambda_atas in the MSbar scheme for a given mass limit
       */
      double getDeltaLambdaAlphatAlphas(int limit, int omitLogs);
      /**
       * 	Returns delta lambda_atas2 in the MSbar scheme for a given mass limit
       * 	@param limit an integer key for a mass limit
       * 	@param omitLogs an integer key to omit all mu terms
       * 	@return delta lambda_atas2 in the MSbar scheme for a given mass limit
       */
      double getDeltaLambdaAlphatAlphas2(int limit, int omitLogs);
      /**
       * 	Returns the Xt^4 terms for Mh^2 at 3-loop level
       * 	@param limit an integer key for a mass limit
       * 	@param omitLogs an integer key to omit all mu terms
       */
      double getXt4Terms(int limit, int omitLogs);
      /**
       * 	Returns the Xt^5 terms for Mh^2 at 3-loop level
       * 	@param limit an integer key for a mass limit
       */
      double getXt5Terms(int limit);
      /**
       * 	Returns the Xt^6 terms for Mh^2 at 3-loop level
       * 	@param limit an integer key for a mass limit
       */
      double getXt6Terms(int limit);
      /**
       * 	Returns the DRbarPrime to MSbar shift of delta lambda 3L at a given xtOrder
       * 	@param limit an integer key for a mass limit
       * 	@param xtOrder an integer key to omit the Xt contributions starting at xtOrder + 1
       * 	@param omitLogs an integer key to omit all log mu terms
       */
      double getDRbarPrimeToMSbarXtTerms(int limit, int xtOrder, int omitLogs);
      /**
       * 	Checks if a threshold correction can be used in the general mass case
       * 	@param exact the threshold correction in the general mass case
       * 	@param shifted the threshold correction in the general mass case with shifted masses
       * 	@param limit the threshold correction in a mass limit
       * 	@returns a bool which is true if the general mass case can be used or false if the limit should be used
       */
      bool isfinite(double exact, double shifted, double limit);
      Parameters p{};	/** The HimalayaInterface struct. */
      double msq2{std::numeric_limits<double>::quiet_NaN()}; /** the average squark mass of the first two generations squared **/
   };
   
}	// himalaya

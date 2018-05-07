#pragma once

#include "Himalaya_interface.hpp"

namespace himalaya{
   
   class ThresholdCalculator{
   public:
      /**
       * 	Constructor
       * 	@param p a HimalayaInterface struct
       * 	@param msq2 the averaged squark mass of the first two generations squared
       * 	@param check a boolean which indicates if the threshold corrections should be tested
       */
      ThresholdCalculator(const Parameters& p, const double msq2, const bool check = false);
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
      Parameters p{};	/** The HimalayaInterface struct. */
      double msq2{};	/** the average squark mass of the first two generations squared **/
   };
   
}	// himalaya
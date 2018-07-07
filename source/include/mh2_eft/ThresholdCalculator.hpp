// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "Himalaya_interface.hpp"
#include <limits>

namespace himalaya{
namespace mh2_eft{

   class ThresholdCalculator{
   public:
      /**
       *         Constructor
       *         @param p_ a HimalayaInterface struct
       *         @param msq2_ the averaged squark mass of the first two generations squared
       *         @param verbose a bool enable the output of the parameter validation. Enabled by default
       *         @param check a boolean which indicates if the threshold corrections should be tested
       */
      ThresholdCalculator(const Parameters& p_, double msq2_ = std::numeric_limits<double>::quiet_NaN(), bool verbose = true, bool check = false);
      /**
       *         Returns a specific threshold corrections for a given mass limit
       *         @param variable an integer key for a threshold correctionn
       *         @param scheme an integer key to set the scheme. Choices are {MSbar, DRbar', DRbar}
       *         @param omitLogs an integer to omit all log mu terms
       *         @return a threshold correction for a given variable in a given scheme for a suitable mass limit
       */
      double getThresholdCorrection(int variable, int scheme, int omitLogs) const;
      /**
       *         Returns the shift needed to convert the 3L threshold correction of lambda to the MSbar scheme
       *         @param xtOrder an integer key to omit the Xt contributions starting at xtOrder + 1
       *         @param omitLogs an integer key to omit all log mu terms
       *         @param omitXtLogs an integer key to omit all Xt^4*Log[mu] and Xt^5*Log[mu] terms
       */
      double getDRbarPrimeToMSbarShift(int xtOrder, int omitLogs, int omitXtLogs = 1) const;
      /**
       *         Returns the DRbarPrime to MSbar shift of delta lambda 3L at a given xtOrder
       *         @param limit an integer key for a mass limit
       *         @param xtOrder an integer key to get a specific Xt order starting at 4
       *         @param omitLogs an integer key to omit all log mu terms
       */
      double getDRbarPrimeToMSbarXtTerms(int limit, int xtOrder, int omitLogs) const;
      /**
       *         Sets the mass limit to check terms
       *         @param limit an integer key for a mass limit
       */
      void setLimit(int limit);
      /**
       *         Get the mass limit determined by ThresholdCalculator
       *         @return The determined mass limit
       */
      int getLimit() const;
      /**
       *         Sets the order of Xt for getDeltaLambdaAlphatAlphas2
       *         @param xtOrder an integer key to truncate delta_lambda_at_as2 at a given Xt order starting at 4
       */
      void setXtOrderOfDeltaLambdaAtAs2(int xtOrder);
   private:
      /**
       *         Returns delta g3_as in the MSbar scheme for a given mass limit
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta g3_as in the MSbar scheme for a given mass limit
       */
      double getDeltaG3Alphas(int omitLogs) const;
      /**
       *         Returns delta yt_as in the MSbar scheme for a given mass limit
       *         @param limit an integer key for a mass limit
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta yt_as in the MSbar scheme for a given mass limit
       */
      double getDeltaYtAlphas(int limit, int omitLogs) const;
      /**
       *         Returns delta yt_as^2 in the MSbar scheme for a given mass limit
       *         @param limit an integer key for a mass limit
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta yt_as^2 in the MSbar scheme for a given mass limit
       */
      double getDeltaYtAlphas2(int limit, int omitLogs) const;
      /**
       *         Returns delta lambda_at in the MSbar scheme for a given mass limit
       *         @param limit an integer key for a mass limit
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta lambda_at in the MSbar scheme for a given mass limit
       */
      double getDeltaLambdaAlphat(int limit, int omitLogs) const;
      /**
       *         Returns delta lambda_atas in the MSbar scheme for a given mass limit
       *         @param limit an integer key for a mass limit
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta lambda_atas in the MSbar scheme for a given mass limit
       */
      double getDeltaLambdaAlphatAlphas(int limit, int omitLogs) const;
      /**
       *         Returns delta lambda_atas2 in the MSbar scheme for a given mass limit
       *         @param limit an integer key for a mass limit
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta lambda_atas2 in the MSbar scheme for a given mass limit
       */
      double getDeltaLambdaAlphatAlphas2(int limit, int omitLogs) const;
      /**
       *         Returns delta_lambda_yb^2_g1^2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yb^2_g1^2
       */
      double getDeltaLambdaYb2G12(int omitLogs) const;
      /**
       *         Returns delta_lambda_g1^4
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_g1^4
       */
      double getDeltaLambdaG14(int omitLogs) const;
      /**
       *         Returns delta_lambda_reg_g1^4
       *         @return delta_lambda_reg_g1^4
       */
      double getDeltaLambdaRegG14() const;
      /**
       *         Returns delta_lambda_chi_g1^4
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_chi_g1^4
       */
      double getDeltaLambdaChiG14(int omitLogs) const;
      /**
       *         Returns delta_lambda_chi_g2^4
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_chi_g2^4
       */
      double getDeltaLambdaChiG24(int omitLogs) const;
      /**
       *         Returns delta_lambda_g2^4
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_g2^4
       */
      double getDeltaLambdaG24(int omitLogs) const;
      /**
       *         Returns delta_lambda_reg_g2^4
       *         @return delta_lambda_reg_g2^4
       */
      double getDeltaLambdaRegG24() const;
      /**
       *         Returns delta_lambda_g1^2_g2^2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_g1^2_g2^2
       */
      double getDeltaLambdaG12G22(int omitLogs) const;
      /**
       *         Returns delta_lambda_reg_g1^2_g2^2
       *         @return delta_lambda_reg_g1^2_g2^2
       */
      double getDeltaLambdaRegG12G22() const;
      /**
       *         Returns delta_lambda_chi_g1^2_g2^2
       *         @return delta_lambda_chi_g1^2_g2^2
       */
      double getDeltaLambdaChiG12G22() const;
      /**
       *         Returns delta_lambda_yb^2_g2^2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yb^2_g2^2
       */
      double getDeltaLambdaYb2G22(int omitLogs) const;
      /**
       *         Returns delta_lambda_yb^4
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yb^4
       */
      double getDeltaLambdaYb4(int omitLogs) const;
      /**
       *         Returns delta_lambda_yt^2_g1^2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yt^2_g1^2
       */
      double getDeltaLambdaYt2G12(int omitLogs) const;
      /**
       *         Returns delta_lambda_yt^2_g2^2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yt^2_g2^2
       */
      double getDeltaLambdaYt2G22(int omitLogs) const;
      /**
       *         Returns delta_lambda_ytau^2_g1^2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_ytau^2_g1^2
       */
      double getDeltaLambdaYtau2G12(int omitLogs) const;
      /**
       *         Returns delta_lambda_ytau^2_g2^2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_ytau^2_g2^2
       */
      double getDeltaLambdaYtau2G22(int omitLogs) const;
      /**
       *         Returns delta_lambda_ytau^4
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_ytau^4
       */
      double getDeltaLambdaYtau4(int omitLogs) const;
      /**
       *         Returns delta_g1_g1
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_g1_g1
       */
      double getDeltaG1G1(int omitLogs) const;
      /**
       *         Returns delta_g2_g2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_g2_g2
       */
      double getDeltaG2G2(int omitLogs) const;
      /**
       *         Returns delta_lambda_v_yt^2
       *         @param limit an integer key for a mass limit
       *         @return delta_lambda_v_yt^2
       */
      double getDeltaVevYt2(int limit) const;
      /**
       *         Returns delta_yt_yt
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_yt_yt
       */
      double getDeltaYtYt(int omitLogs) const;
      /**
       *         Returns delta_ytau_ytau
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_ytau_ytau
       */
      double getDeltaYtauYtau(int omitLogs) const;
      /**
       *         Returns delta_yt_yb
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_yt_yb
       */
      double getDeltaYtYb(int omitLogs) const;
      /**
       *         Returns delta_lambda_yb4_g32
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yb4_g32
       */
      double getDeltaLambdaYb4G32(int omitLogs) const;
      /**
       *         Returns delta_lambda_yb6
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yb6
       */
      double getDeltaLambdaYb6(int omitLogs) const;
      /**
       *         Returns delta_lambda_yt6
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yt6
       */
      double getDeltaLambdaYt6(int omitLogs) const;
      /**
       *         Returns delta_lambda_ytau6
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_ytau6
       */
      double getDeltaLambdaYtau6(int omitLogs) const;
      /**
       *         Returns delta_lambda_yt2_yb4
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yt2_yb4
       */
      double getDeltaLambdaYt2Yb4(int omitLogs) const;
      /**
       *         Returns delta_lambda_yt4_yb2
       *         @param omitLogs an integer key to omit all mu terms
       *         @return delta_lambda_yt4_yb2
       */
      double getDeltaLambdaYt4Yb2(int omitLogs) const;
      Parameters p{}; ///< The HimalayaInterface struct
      double msq2{std::numeric_limits<double>::quiet_NaN()}; ///< the average squark mass of the first two generations squared
      int xtOrderLambdaAtAs2 = 6; ///< A flag to truncate the Xt order of delta_lambda_at_as2 at a given value starting at 4
   };

}        // mh2_eft
}        // himalaya

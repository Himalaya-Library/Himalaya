// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

/**
 * @file EFTFlags.hpp
 *
 * @brief Enum definitions for the EFT calculation
 */

namespace himalaya{
namespace mh2_eft{


/// Renormalization scheme numbers
namespace RenSchemes {
   /// renormalization scheme falgs
   enum RenSchems : int {
      FIRST = 0,                    ///< [0]
      H3m = FIRST,                  ///< [0] H3m
      DRBARPRIME,                   ///< [1] DRbar'
      H3mMDRBAR,                    ///< [2] MDRbar with H3m renormalization
      MDRBARPRIME,                  ///< [3] MDRbar with DRbar' renormalization
      NUMBER_OF_REN_SCHEMES,        ///< [4] number of actual renormalization schemes
      TEST = NUMBER_OF_REN_SCHEMES, ///< [4] A non-physical scheme to test implemented expressions
      NUMBER_OF_ALL_REN_SCHEMES     ///< [5] number of all renormalization schemes
   };
} // namespace RenSchemes


namespace ThresholdVariables {
   /// Threshold variable flags
   enum ThresholdVariables : int {
      G3_AS,              ///< g3_as threshold correction
      YT_AS,              ///< yt_as threshold correction
      YT_AS2,             ///< yt_as^2 threshold correction
      LAMBDA_AT,          ///< lambga_at threshold correction
      LAMBDA_AT_AS,       ///< lambda_atas threshold correction
      LAMBDA_AT_AS2,      ///< lambda_atas2 threshold correction
      LAMBDA_YB2_G12,     ///< lambda_yb2_g12 threshold correction
      LAMBDA_G14,         ///< lambda_g14 threshold correction
      LAMBDA_REG_G14,     ///< lambda_reg_g14 threshold correction
      LAMBDA_CHI_G14,     ///< lambda_chi_g14 threshold correction
      LAMBDA_CHI_G24,     ///< lambda_chi_g24 threshold correction
      LAMBDA_G24,         ///< lambda_g24 threshold correction
      LAMBDA_REG_G24,     ///< lambda_reg_g24 threshold correction
      LAMBDA_G12_G22,     ///< lambda_g12_g22 threshold correction
      LAMBDA_REG_G12_G22, ///< lambda_reg_g12_g22 threshold correction
      LAMBDA_CHI_G12_G22, ///< lambda_chi_g12_g22 threshold correction
      LAMBDA_YB2_G22,     ///< lambda_yb2_g22 threshold correction
      LAMBDA_YB4,         ///< lambda_yb4 threshold correction
      LAMBDA_YT2_G12,     ///< lambda_yt2_g12 threshold correction
      LAMBDA_YT2_G22,     ///< lambda_yt2_g22 threshold correction
      LAMBDA_YTAU2_G12,   ///< lambda_ytau2_g12 threshold correction
      LAMBDA_YTAU2_G22,   ///< lambda_ytau2_g22 threshold correction
      LAMBDA_YTAU4,       ///< lambda_ytau4 threshold correction
      G1_G1,              ///< g1_g1 threshold correction
      G2_G2,              ///< g2_g2 threshold correction
      VEV_YT2,            ///< vev_yt2 threshold correction
      YT_YB,              ///< yt_yb threshold correction
      YT_YT,              ///< yt_yt threshold correction
      YTAU_YTAU,          ///< ytau_ytau threshold correction
      LAMBDA_YB4_G32,     ///< yb4_g32 threshold correction
      LAMBDA_YB6,         ///< yb6s threshold correction
      LAMBDA_YT6,         ///< yt6 threshold correction
      LAMBDA_YTAU6,       ///< ytau6 threshold correction
      LAMBDA_YT2_YB4,     ///< yt2_yb4 threshold correction
      LAMBDA_YT4_YB2,     ///< yt4_yb2 threshold correction
      VEV_YB2,            ///< vev_yb2 threshold correction
      VEV_YTAU2,          ///< vev_ytau2 threshold correction
      VEV_G12,            ///< ytau_yb threshold correction
      VEV_G22,            ///< ytau_yb threshold correction
      YTAU_YB,            ///< ytau_yb threshold correction
      LAMBDA_YTAU4_YB2,   ///< ytau4_yb2 threshold correction
      LAMBDA_YTAU2_YB4,   ///< ytau2_yb4 threshold correction
      YB_YT,
      YB_AS,
      YB_YB,
      NUMBER_OF_THRESHOLD_VARIALES
   };
} // namespace ThresholdVariables


namespace EFTOrders {
   /// Order flags
   enum EFTOrders : int{
      G12G22,
      G14,
      G24,
      G12YB2,
      G22YB2,
      YB4,
      G12YTAU2,
      G22YTAU2,
      YTAU4,
      G12YT2,
      G22YT2,
      YT4,
      G32YB4,
      YB6,
      YT6,
      YTAU2YB4,
      YTAU2YT4,//??
      YTAU6,
      YT2YB4,
      G32YT4,
      YB2YT4,
      YTAU4YB2,
      ONLY_AT_AS,
      NUMBER_OF_EFT_ORDERS
   };
} // namespace EFTOrders


/// limits
enum class Limits : int {
   GENERAL,          ///< general mass case
   MQ3_EQ_MU3,       ///< mQ3 = mU3
   MQ3_EQ_M3,        ///< mQ3 = m3
   MU3_EQ_M3,        ///< mU3 = m3
   MQ3_EQ_MU3_EQ_M3, ///< mQ3 = mU3 = m3
   DEGENERATE,       ///< mQ3 = mU3 = m3 = msq
   NUMBER_OF_LIMITS
};

}        // mh2_eft
}        // himalaya

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
      H3m           = 0,        ///< H3m
      DRBARPRIME    = 1,        ///< DRbar'
      H3mMDRBAR     = 2,        ///< MDRbar with H3m renormalization
      MDRBARPRIME   = 3,        ///< MDRbar with DRbar' renormalization
      TEST          = 4,        ///< A non-physical scheme to test implemented expressions
      NUMBER_OF_REN_SCHEMES
   };
} // namespace RenSchemes

namespace ThresholdVariables {

   /// Threshold variable flags
   enum ThresholdVariables : int{
      G3_AS              = 0,       ///< g3_as threshold correction
      YT_AS              = 1,       ///< yt_as threshold correction
      YT_AS2             = 2,       ///< yt_as^2 threshold correction
      LAMBDA_AT          = 3,       ///< lambga_at threshold correction
      LAMBDA_AT_AS       = 4,       ///< lambda_atas threshold correction
      LAMBDA_AT_AS2      = 5,       ///< lambda_atas2 threshold correction
      LAMBDA_YB2_G12     = 6,       ///< lambda_yb2_g12 threshold correction
      LAMBDA_G14         = 7,       ///< lambda_g14 threshold correction
      LAMBDA_REG_G14     = 8,       ///< lambda_reg_g14 threshold correction
      LAMBDA_CHI_G14     = 9,       ///< lambda_chi_g14 threshold correction
      LAMBDA_CHI_G24     = 10,      ///< lambda_chi_g24 threshold correction
      LAMBDA_G24         = 11,      ///< lambda_g24 threshold correction
      LAMBDA_REG_G24     = 12,      ///< lambda_reg_g24 threshold correction
      LAMBDA_G12_G22     = 13,      ///< lambda_g12_g22 threshold correction
      LAMBDA_REG_G12_G22 = 14,      ///< lambda_reg_g12_g22 threshold correction
      LAMBDA_CHI_G12_G22 = 15,      ///< lambda_chi_g12_g22 threshold correction
      LAMBDA_YB2_G22     = 16,      ///< lambda_yb2_g22 threshold correction
      LAMBDA_YB4         = 17,      ///< lambda_yb4 threshold correction
      LAMBDA_YT2_G12     = 18,      ///< lambda_yt2_g12 threshold correction
      LAMBDA_YT2_G22     = 19,      ///< lambda_yt2_g22 threshold correction
      LAMBDA_YTAU2_G12   = 20,      ///< lambda_ytau2_g12 threshold correction
      LAMBDA_YTAU2_G22   = 21,      ///< lambda_ytau2_g22 threshold correction
      LAMBDA_YTAU4       = 22,      ///< lambda_ytau4 threshold correction
      G1_G1              = 23,      ///< g1_g1 threshold correction
      G2_G2              = 24,      ///< g2_g2 threshold correction
      VEV_YT2            = 25,      ///< vev_yt2 threshold correction
      YT_YB              = 26,      ///< yt_yb threshold correction
      YT_YT              = 27,      ///< yt_yt threshold correction
      YTAU_YTAU          = 28,      ///< ytau_ytau threshold correction
      LAMBDA_YB4_G32     = 29,      ///< yb4_g32 threshold correction
      LAMBDA_YB6         = 30,      ///< yb6s threshold correction
      LAMBDA_YT6         = 31,      ///< yt6 threshold correction
      LAMBDA_YTAU6       = 32,      ///< ytau6 threshold correction
      LAMBDA_YT2_YB4     = 33,      ///< yt2_yb4 threshold correction
      LAMBDA_YT4_YB2     = 34,      ///< yt4_yb2 threshold correction
      VEV_YB2            = 35,      ///< vev_yb2 threshold correction
      VEV_YTAU2          = 36,      ///< vev_ytau2 threshold correction
      VEV_G12            = 37,      ///< ytau_yb threshold correction
      VEV_G22            = 38,      ///< ytau_yb threshold correction
      YTAU_YB            = 39,      ///< ytau_yb threshold correction
      LAMBDA_YTAU4_YB2   = 40,      ///< ytau4_yb2 threshold correction
      LAMBDA_YTAU2_YB4   = 41,      ///< ytau2_yb4 threshold correction
      YB_YT              = 42,
      YB_AS              = 43,
      YB_YB              = 44,
      NUMBER_OF_THRESHOLD_VARIALES
   };
} // namespace ThresholdVariables

namespace EFTOrders {

   /// Order flags
   enum EFTOrders : int{
      FIRST    = 0,
      G12G22   = FIRST,
      G14      = 1,
      G24      = 2,
      G12YB2   = 3,
      G22YB2   = 4,
      YB4      = 5,
      G12YTAU2 = 6,
      G22YTAU2 = 7,
      YTAU4    = 8,
      G12YT2   = 9,
      G22YT2   = 10,
      YT4      = 11,
      G32YB4   = 12,
      YB6      = 13,
      YT6      = 14,
      YTAU2YB4 = 15,
      YTAU2YT4 = 16,//??
      YTAU6    = 17,
      YT2YB4   = 18,
      G32YT4   = 19,
      YB2YT4   = 20,
      YTAU4YB2 = 21,
      ONLY_AT_AS = 22,
      NUMBER_OF_EFT_ORDERS
   };
} // namespace EFTOrders

/// Limit numbers
namespace Limits {

   /// limit flags
   enum Limits : int{
      GENERAL                 = 0,      ///< general mass case
      MQ3_EQ_MU3              = 1,      ///< mQ3 = mU3
      MQ3_EQ_M3               = 2,      ///< mQ3 = m3
      MU3_EQ_M3               = 3,      ///< mU3 = m3
      MQ3_EQ_MU3_EQ_M3        = 4,      ///< mQ3 = mU3 = m3
      DEGENERATE              = 5,      ///< mQ3 = mU3 = m3 = msq
      NUMBER_OF_LIMITS
   };
} // namespace Limits

}        // mh2_eft
}        // himalaya

#pragma once

namespace himalaya {
namespace Hierarchies {

/// Hierarchy numbers
enum Hierarchies : int {
   FIRST   =  0,
   h3      =  FIRST,        ///< The key to hierarchy h3
   h32q2g  =  1,            ///< The key to hierarchy h32q2g
   h3q22g  =  2,            ///< The key to hierarchy h3q22g
   h4      =  3,            ///< The key to hierarchy h4
   h5      =  4,            ///< The key to hierarchy h5
   h5g1    =  5,            ///< The key to hierarchy h5g1
   h6      =  6,            ///< The key to hierarchy h6
   h6b     =  7,            ///< The key to hierarchy h6b
   h6b2qg2 =  8,            ///< The key to hierarchy h6b2qg2
   h6bq22g =  9,            ///< The key to hierarchy h6bq22g
   h6bq2g2 = 10,            ///< The key to hierarchy h6bq2g2
   h6g2    = 11,            ///< The key to hierarchy h6g2
   h9      = 12,            ///< The key to hierarchy h9
   h9q2    = 13,            ///< The key to hierarchy h9q2
   NUMBER_OF_HIERARCHIES    ///< Number of hierarchies
};

} // namespace Hierarchies

namespace ExpansionDepth {

/// expansion depth flags
enum ExpansionDepth : int {
   FIRST     = 14,
   xx        = FIRST,   ///< This flag can truncate the two loop expansion at the three loop expansion depth
   xxMst     = 15,      ///< This flag can truncate the expansion depth of the stop/sbottom masses by one order
   xxDmglst1 = 16,      ///< This flag can truncate the expansion depth of the difference of stop/sbottom 1 mass and the gluino mass by one order
   xxDmsqst1 = 17,      ///< This flag can truncate the expansion depth of the difference of the stop/sbottom 1 mass and the average squark mass by one order
   xxDmst12  = 18,      ///< This flag can truncate the expansion depth of the difference of the stop/sbottom masses by one order
   xxAt      = 19,      ///< This flag can truncate the expansion depth of At/Ab by one order
   xxlmMsusy = 20,      ///< This flag can truncate the expansion depth of log(Msusy) by one order
   xxMsq     = 21,      ///< This flag can truncate the expansion depth of the average squark mass by one order
   xxMsusy   = 22,      ///< This flag can truncate the expansion depth of the average SUSY mass by one order
   xxDmglst2 = 23,      ///< This flag can truncate the expansion depth of the difference of the stop/sbottom 2 mass and the gluino mass by one order
   xxDmsqst2 = 24,      ///< This flag can truncate the expansion depth of the difference of the average squark mass and the stop/sbottom 2 mass by one order
   xxMgl     = 25,      ///< This flag can truncate the expansion depth of the gluino mass by one order
   NUMBER_OF_EXPANSIONS ///< Number of expansions
};

} // namespace ExpansionDepth

/// Mass schemes
namespace MassSchemes{
   
   /// mass scheme flags
   enum MassSchemes : int{
      FIRST        = 0,
      DEFAULT      = FIRST,
      MASSEIGEN    = 1,
      SOFTBREAKING = 2,
      NUMBER_OF_SCHEMES
   };
}

/// Limit numbers
namespace Limits{
   
   /// limit flags
   enum Limits : int{
      FIRST                   = 0,
      GENERAL                 = FIRST,  ///< Flag for general mass case
      MQ3_EQ_MU3              = 1,      ///< Flag for mQ3 = mU3
      MQ3_EQ_M3               = 2,      ///< Flag for mQ3 = m3
      MU3_EQ_M3               = 3,      ///< Flag for mU3 = m3
      MQ3_EQ_MU3_EQ_M3        = 4,      ///< Flag for mQ3 = mU3 = m3
      DEGENERATE              = 5,	///< Flag for mQ3 = mU3 = m3 = msq
      NUMBER_OF_LIMITS                  ///< Number of limits
   };
} // namespace Limits

/// Renormalization scheme numbers
namespace RenSchemes{
   
   /// renormalization scheme falgs
   enum RenSchems : int {
      FIRST         = 0,
      H3m           = FIRST,	///< H3m flag
      DRBARPRIME    = 1,	///< DRbar' flag
      H3mMDRBAR     = 2,	///< MDRbar flag with H3m renormalization
      MDRBARPRIME   = 3,	///< MDRbar flag with DRbar' renormalization
      TEST          = 4,	///< A non-physical scheme to test implemented expressions
      NUMBER_OF_REN_SCHEMES	///< Number of renormalization schemes
   };
} // namespace RenSchemes

namespace ThresholdVariables{
   
   /// Threshold variable flags
   enum ThresholdVariables : int{
      FIRST         = 0,
      G3_AS         = FIRST,        ///< Flag for g3_as threshold correction
      YT_AS         = 1,            ///< Flag for yt_as threshold correction
      YT_AS2        = 2,            ///< Flag for yt_as^2 threshold correction
      LAMBDA_AT     = 3,            ///< Flag for lambga_at threshold correction
      LAMBDA_AT_AS  = 4,            ///< Flag for lambda_atas threshold correction
      LAMBDA_AT_AS2 = 5,            ///< Flag for lambda_atas2 threshold correction
      NUMBER_OF_THRESHOLD_VARIALES  ///< Number of threshold variables
   };
} // namespace ThresholdVariables

namespace EFTOrders{
   
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
      ALL      = 21,
      NUMBER_OF_EFT_ORDERS
   };
}

} // namespace himalaya

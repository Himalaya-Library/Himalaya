// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

namespace himalaya {
namespace hierarchies {

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

} // namespace hierarchies
} // namespace himalaya

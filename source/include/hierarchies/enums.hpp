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
   h3      =  FIRST,
   h32q2g  =  1,
   h3q22g  =  2,
   h4      =  3,
   h5      =  4,
   h5g1    =  5,
   h6      =  6,
   h6b     =  7,
   h6b2qg2 =  8,
   h6bq22g =  9,
   h6bq2g2 = 10,
   h6g2    = 11,
   h9      = 12,
   h9q2    = 13,
   NUMBER_OF_HIERARCHIES
};

} // namespace Hierarchies

namespace ExpansionDepth {

/// expansion depth flags
enum ExpansionDepth : int {
   FIRST     =  0,
   threeLoop =  FIRST,   ///< truncate the two loop expansion at the three loop expansion depth
   Mst       =  1,       ///< truncate the expansion depth in the stop/sbottom masses by one order
   Dmglst1   =  2,       ///< truncate the expansion depth in the difference of stop/sbottom 1 mass and the gluino mass by one order
   Dmsqst1   =  3,       ///< truncate the expansion depth in the difference of the stop/sbottom 1 mass and the average squark mass by one order
   Dmst12    =  4,       ///< truncate the expansion depth in the difference of the stop/sbottom masses by one order
   At        =  5,       ///< truncate the expansion depth in At/Ab by one order
   lmMsusy   =  6,       ///< truncate the expansion depth in log(Msusy) by one order
   Msq       =  7,       ///< truncate the expansion depth in the average squark mass by one order
   Msusy     =  8,       ///< truncate the expansion depth in the average SUSY mass by one order
   Dmglst2   =  9,       ///< truncate the expansion depth in the difference of the stop/sbottom 2 mass and the gluino mass by one order
   Dmsqst2   = 10,       ///< truncate the expansion depth in the difference of the average squark mass and the stop/sbottom 2 mass by one order
   Mgl       = 11,       ///< truncate the expansion depth in the gluino mass by one order
   NUMBER_OF_EXPANSIONS
};

} // namespace ExpansionDepth

/// Mass schemes
namespace MassSchemes {

   /// mass scheme flags
   enum MassSchemes : int {
      FIRST        = 0,
      DEFAULT      = FIRST,
      MASSEIGEN    = 1,
      SOFTBREAKING = 2,
      NUMBER_OF_SCHEMES
   };
}

} // namespace hierarchies
} // namespace himalaya

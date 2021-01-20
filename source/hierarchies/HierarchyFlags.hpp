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
   threeLoop,       ///< truncate the two loop expansion at the three loop expansion depth
   Mst,             ///< truncate the expansion depth in the stop/sbottom masses by one order
   Dmglst1,         ///< truncate the expansion depth in the difference of stop/sbottom 1 mass and the gluino mass by one order
   Dmsqst1,         ///< truncate the expansion depth in the difference of the stop/sbottom 1 mass and the average squark mass by one order
   Dmst12,          ///< truncate the expansion depth in the difference of the stop/sbottom masses by one order
   At,              ///< truncate the expansion depth in At/Ab by one order
   lmMsusy,         ///< truncate the expansion depth in log(Msusy) by one order
   Msq,             ///< truncate the expansion depth in the average squark mass by one order
   Msusy,           ///< truncate the expansion depth in the average SUSY mass by one order
   Dmglst2,         ///< truncate the expansion depth in the difference of the stop/sbottom 2 mass and the gluino mass by one order
   Dmsqst2,         ///< truncate the expansion depth in the difference of the average squark mass and the stop/sbottom 2 mass by one order
   Mgl,             ///< truncate the expansion depth in the gluino mass by one order
   NUMBER_OF_EXPANSIONS
};

} // namespace ExpansionDepth

} // namespace hierarchies
} // namespace himalaya

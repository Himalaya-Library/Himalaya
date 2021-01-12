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
   h3      =  FIRST,        ///< hierarchy h3
   h32q2g  =  1,            ///< hierarchy h32q2g
   h3q22g  =  2,            ///< hierarchy h3q22g
   h4      =  3,            ///< hierarchy h4
   h5      =  4,            ///< hierarchy h5
   h5g1    =  5,            ///< hierarchy h5g1
   h6      =  6,            ///< hierarchy h6
   h6b     =  7,            ///< hierarchy h6b
   h6b2qg2 =  8,            ///< hierarchy h6b2qg2
   h6bq22g =  9,            ///< hierarchy h6bq22g
   h6bq2g2 = 10,            ///< hierarchy h6bq2g2
   h6g2    = 11,            ///< hierarchy h6g2
   h9      = 12,            ///< hierarchy h9
   h9q2    = 13,            ///< hierarchy h9q2
   NUMBER_OF_HIERARCHIES    ///< Number of hierarchies
};

} // namespace Hierarchies

namespace ExpansionDepth {

/// expansion depth flags
enum ExpansionDepth : int {
   FIRST     =  0,
   xx        =  FIRST,   ///< truncate the two loop expansion at the three loop expansion depth
   xxMst     =  1,       ///< truncate the expansion depth in the stop/sbottom masses by one order
   xxDmglst1 =  2,       ///< truncate the expansion depth in the difference of stop/sbottom 1 mass and the gluino mass by one order
   xxDmsqst1 =  3,       ///< truncate the expansion depth in the difference of the stop/sbottom 1 mass and the average squark mass by one order
   xxDmst12  =  4,       ///< truncate the expansion depth in the difference of the stop/sbottom masses by one order
   xxAt      =  5,       ///< truncate the expansion depth in At/Ab by one order
   xxlmMsusy =  6,       ///< truncate the expansion depth in log(Msusy) by one order
   xxMsq     =  7,       ///< truncate the expansion depth in the average squark mass by one order
   xxMsusy   =  8,       ///< truncate the expansion depth in the average SUSY mass by one order
   xxDmglst2 =  9,       ///< truncate the expansion depth in the difference of the stop/sbottom 2 mass and the gluino mass by one order
   xxDmsqst2 = 10,       ///< truncate the expansion depth in the difference of the average squark mass and the stop/sbottom 2 mass by one order
   xxMgl     = 11,       ///< truncate the expansion depth in the gluino mass by one order
   NUMBER_OF_EXPANSIONS
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

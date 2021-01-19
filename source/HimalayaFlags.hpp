// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

/**
 * @file HimalayaFlags.hpp
 *
 * @brief enum definitions
 */

namespace himalaya {

namespace CouplingOrders {
   /// Coupling orders for calculation
   enum CouplingOrders : int {
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
} // namespace CouplingOrders

} // namespace himalaya

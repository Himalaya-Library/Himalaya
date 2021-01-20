// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

/**
 * @file pv.hpp
 *
 * @brief Declaration of real Passarino-Veltman loop functions with
 * squared arguments.
 */

#pragma once

namespace himalaya {
namespace mh2_fo {

/// A0 Passarino-Veltman function
double a0(double m2, double q2) noexcept;
/// B0 Passarino-Veltman function
double b0(double p2, double m12, double m22, double q2) noexcept;
/// derivative of B0 Passarino-Veltman function w.r.t. p^2, for p^2 = 0
double d1_b0(double m12, double m22) noexcept;

} // namespace mh2_fo
} // namespace himalaya

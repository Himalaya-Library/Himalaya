// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

/**
 * @file pv.hpp
 *
 * @brief Real Passarino-Veltman loop functions with squared
 * arguments.
 */

#pragma once

namespace himalaya {
namespace mh2_fo {

double a0(double m2, double q2) noexcept;
double b0(double p2, double m12, double m22, double q2) noexcept;
double d1_b0(double m1, double m2) noexcept;

} // namespace mh2_fo
} // namespace himalaya

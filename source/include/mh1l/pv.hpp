// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

/**
 * @file pv.hpp
 *
 * @brief Real Passarino-Veltman loop functions with squared
 * arguments.
 */

#pragma once

namespace himalaya {
namespace mh1l {

double a0(double m2, double q2) noexcept;
double b0(double p2, double m12, double m22, double q2) noexcept;

} // namespace mh1l
} // namespace himalaya

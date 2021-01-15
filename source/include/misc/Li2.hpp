// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#pragma once

#include <complex>

/**
 * @file Li2.hpp
 * @brief Declaration of the dilogarithm function
 * @note The file has been taken from GM2Calc.
 */

#ifdef __GNUC__
#  define ATTR(x) __attribute__ ((x))
#else
#  define ATTR(x)
#endif

namespace himalaya {

/// real dilogarithm
double dilog(double) noexcept ATTR(const);

/// complex dilogarithm
std::complex<double> dilog(const std::complex<double>&) noexcept ATTR(const);

/// Clausen function Cl_2(x)
double clausen_2(double) noexcept ATTR(const);

} // namespace himalaya

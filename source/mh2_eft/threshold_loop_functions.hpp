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

#pragma once

/**
 * @file threshold_loop_functions.hpp
 * @brief Declaration of the loop functions from [arXiv:1407.4081]
 * @note The file has been taken from FlexibleSUSY.
 */

namespace himalaya {
namespace mh2_eft {
namespace threshold_loop_functions {

#ifdef __GNUC__
#  define ATTR(x) __attribute__ ((x))
#else
#  define ATTR(x)
#endif

#define TCFATTR noexcept ATTR(const)

double F1(double) TCFATTR;
double F2(double) TCFATTR;
double F3(double) TCFATTR;
double F4(double) TCFATTR;
double F5(double) TCFATTR;
double F6(double) TCFATTR;
double F7(double) TCFATTR;
double F8(double, double) TCFATTR;
double F9(double, double) TCFATTR;

double f(double) TCFATTR;
double g(double) TCFATTR;

double f1(double) TCFATTR;
double f2(double) TCFATTR;
double f3(double) TCFATTR;
double f4(double) TCFATTR;
double f5(double, double) TCFATTR;
double f6(double, double) TCFATTR;
double f7(double, double) TCFATTR;
double f8(double, double) TCFATTR;

// 2-loop threshold function fth[y] from MhEFT-1.1
double fth1(double) TCFATTR;
double fth2(double) TCFATTR;
double fth3(double) TCFATTR;

/// \f$I_{abc}(a,b,c)\f$ (arguments are interpreted as unsquared)
double Iabc(double, double, double) TCFATTR;

/// \f$Delta_{xyz}(x,y,z)\f$ (arguments are interpreted as squared masses)
double delta_xyz(double, double, double) TCFATTR;

/// \f$phi_{xyz}(x,y,z)\f$ (arguments are interpreted as squared masses)
double phi_xyz(double, double, double) TCFATTR;

} // namespace threshold_loop_functions
} // namespace mh2_eft
} // namespace himalaya

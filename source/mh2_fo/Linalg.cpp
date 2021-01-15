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

#include "Linalg.hpp"
#include <cmath>

/**
 * @file Linalg.cpp
 * @brief Implementation of the Passarino-Veltman functions
 * @note This file has been taken from FlexibleSUSY.
 */

namespace himalaya {

namespace {

double pow2(double x) noexcept { return x*x; }

} // anonymous namespace

/**
 * Diagonalizes a mass matrix perturbatively.
 *
 * @param m0 tree-level contribution
 * @param m1 1-loop contribution
 * @param m2 2-loop contribution
 * @param m3 3-loop contribution
 *
 * @return perturbatively calculated mass eigenvalues
 */
std::tuple<V2,V2,V2,V2> fs_diagonalize_hermitian_perturbatively(
   const RM22& m0, const RM22& m1, const RM22& m2, const RM22& m3)
{
   // tree-level
   const auto a11 = m0(0,0), a12 = m0(0,1), a22 = m0(1,1);
   const auto c1 = a11 + a22;
   const auto c2 = std::sqrt(pow2(a11) + 4*pow2(a12) - 2*a11*a22 + pow2(a22));

   V2 mh2_0L;
   mh2_0L << 0.5*(c1 - c2), 0.5*(c1 + c2);

   // 1-loop
   const auto b11 = m1(0,0), b12 = m1(0,1), b22 = m1(1,1);
   const auto c3 = b11 + b22;
   const auto c4 = (a11*b11 - a22*b11 + 4*a12*b12 - a11*b22 + a22*b22)/c2;

   V2 mh2_1L;
   mh2_1L << 0.5*(c3 - c4), 0.5*(c3 + c4);

   // 2-loop
   const auto d11 = m2(0,0), d12 = m2(0,1), d22 = m2(1,1);
   const auto c5 = d11 + d22;
   const auto c6 = 0.5*(pow2(b11) + 4*pow2(b12) - 2*b11*b22 + pow2(b22)
                        - pow2(c4) + 2*a11*d11 - 2*a22*d11 + 8*a12*d12
                        - 2*a11*d22 + 2*a22*d22)/c2;

   V2 mh2_2L;
   mh2_2L << 0.5*(c5 - c6), 0.5*(c5 + c6);

   // 3-loop
   const auto e11 = m3(0,0), e12 = m3(0,1), e22 = m3(1,1);
   const auto c7 = e11 + e22;
   const auto c8 = 0.5/c2*(
      2*(4*b12*d12 + b11*(d11 - d22) + b22*(-d11 + d22) + a11*e11 - a22*e11
         + 4*a12*e12 - a11*e22 + a22*e22)
      + (c4*(2*b11*b22 - 2*a11*d11 + 2*a22*d11 - 8*a12*d12 + 2*a11*d22
             - 2*a22*d22 - pow2(b11) - 4*pow2(b12) - pow2(b22) + pow2(c4)))/c2);

   V2 mh2_3L;
   mh2_3L << 0.5*(c7 - c8), 0.5*(c7 + c8);

   return std::make_tuple(mh2_0L, mh2_1L, mh2_2L, mh2_3L);
}

} // namespace himalaya

// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

namespace himalaya {
namespace {

const double Pi        = 3.14159265358979324;      ///< Pi
const double z2        = 1.64493406684822644;      ///< Zeta[2] = Pi^2/6
const double z3        = 1.20205690315959429;      ///< Zeta[3]
const double z4        = 1.08232323371113819;      ///< Zeta[4] = Pi^4/90
const double sqrt2     = 1.41421356237309505;      ///< Sqrt[2]
const double sqrt3     = 1.73205080756887729;      ///< Sqrt[3]
const double sqrt15    = 3.87298334620741689;      ///< Sqrt[15]
const double sqrt35    = 0.774596669241483377;     ///< Sqrt[3/5]
const double inv_sqrt2 = 0.707106781186547524;     ///< 1/Sqrt[2]
const double log2      = 0.693147180559945309;     ///< Log[2]
const double oneLoop   = 6.332573977646110963e-03; ///< 1/(4 Pi)^2
const double twoLoop   = 4.010149318236068752e-05; ///< 1/(4 Pi)^4
const double threeLoop = 2.539456721913701978e-07; ///< 1/(4 Pi)^6
const double B4        = -1.7628000870737709;      ///< -4 Zeta[2] Log[2]^2 + 2/3 Log[2]^4 - 13/2 Zeta[4] + 16 PolyLog[4,1/2]
const double D3        = -3.0270094939876520;      ///< 6 Zeta[3] - 15/4 Zeta[4] - 6 Im[PolyLog[2, Exp[I Pi/3]]]^2
const double DN        =  1.1202483970392421;      ///< 6 Zeta[3] - 4 Zeta[2] Log[2]^2 + 2/3 Log[2]^4 - 21/2 Zeta[4] + 16 PolyLog[4,1/2]
const double OepS2     =  7.8517428364255312;      ///< -763/32 - 9 Pi Sqrt[3] Log[3]^2/16 - 35 Pi^3 Sqrt[3]/48 + 195/16 Zeta[2] - 15/4 Zeta[3] + 57/16 Zeta[4] + 45 Sqrt[3]/2 Im[PolyLog[2, Exp[I Pi / 3]]] - 27 Sqrt[3] Im[PolyLog[3, Exp[-I Pi/6]/Sqrt[3]]]
const double S2        = 0.26043413763216210;      ///< 4 Im[PolyLog[2, Exp[I Pi/3]]]/(9 Sqrt[3])
const double T1ep      = -24.208928021203593;      ///< -45/2 - Pi Sqrt[3] Log[3]^2/8 - 35 Pi^3 Sqrt[3]/216 - 9/2 Zeta[2] + Zeta[3] + 6 Sqrt[3] Im[PolyLog[2, Exp[I Pi/3]]] - 6 Sqrt[3] Im[PolyLog[3, Exp[-I Pi/6]/Sqrt[3]]]

} // anonymous namespace
} // namespace himalaya

// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "power.hpp"
#include <complex>

namespace himalaya {
namespace {

const double Pi        = 3.14159265358979324;      ///< Pi
const double z2        = 1.64493406684822644;      ///< Zeta[2] = Pi^2/6
const double z3        = 1.20205690315959429;      ///< Zeta[3]
const double z4        = 1.08232323371113819;      ///< Zeta[4] = Pi^4/90
const double sqrt2     = 1.41421356237309505;      ///< Sqrt[2]
const double sqrt3     = 1.73205080756887729;      ///< Sqrt[3]
const double log2      = 0.693147180559945309;     ///< Log[2]
const double log3      = 1.09861228866810969;      ///< Log[3]
const double oneLoop   = 6.332573977646110963e-03; ///< 1/(4 Pi)^2
const double twoLoop   = 4.010149318236068752e-05; ///< 1/(4 Pi)^4
const double threeLoop = 2.539456721913701978e-07; ///< 1/(4 Pi)^6

// polylogs
const double pl412 = 0.51747906167389934317668576113647; ///< PolyLog[4,1/2]
const std::complex<double> pl2expPi3(0.27415567780803773941206919444, 1.014941606409653625021202554275); ///< PolyLog[2, Exp[I Pi / 3]]
const std::complex<double> pl3expPi6sqrt3(0.51928806536375962552715984277228, - 0.33358157526196370641686908633664); ///< PolyLog[3, Exp[- I Pi / 6] / Sqrt[3]]

// polylog functions
const double B4 = -4 * z2 * pow2(log2) + 2 / 3.* pow4(log2) - 13 / 2. * z4 + 16. * pl412;
const double D3 = 6 * z3 - 15 / 4. * z4 - 6. * pow2(std::imag(pl2expPi3));
const double DN = 6 * z3 - 4 * z2 * pow2(log2) + 2 / 3. * pow4(log2) - 21 / 2. * z4 + 16. * pl412;
const double OepS2 = - 763 / 32. - (9 * Pi * sqrt3 * pow2(log3)) / 16. - (35 * pow3(Pi) * sqrt3) / 48.
   + 195 / 16. * z2 - 15 / 4. * z3 + 57 / 16. * z4 + 45 * sqrt3 / 2. * std::imag(pl2expPi3)
   - 27 * sqrt3 * std::imag(pl3expPi6sqrt3);
const double S2 = 4 * std::imag(pl2expPi3) / (9. * sqrt3);
const double T1ep = - 45 / 2. - (Pi * sqrt3 * pow2(log3)) / 8. - (35 * pow3(Pi) * sqrt3) / 216. - 9 / 2. * z2 + z3 
   + 6. * sqrt3 * std::imag(pl2expPi3) - 6. * sqrt3 * std::imag(pl3expPi6sqrt3);

} // anonymous namespace
} // namespace himalaya

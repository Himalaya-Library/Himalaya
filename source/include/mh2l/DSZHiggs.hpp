// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <Eigen/Core>

/**
 * @file DSZHiggs.hpp
 * @brief function declarations for 2-loop MSSM Higgs contributions
 *
 * Notation:
 *
 * mt2    : squared DR' top mass in the MSSM
 * mb2    : squared DR' bottom mass in the MSSM
 * mtau2  : squared DR' tau mass in the MSSM
 * mg     : DR' gluino mass in the MSSM
 * mA2    : squared DR' CP-odd Higgs mass in the MSSM
 * mst12  : squared DR' lightest stop mass
 * mst22  : squared DR' heaviest stop mass
 * msb12  : squared DR' lightest sbottom mass
 * msb22  : squared DR' heaviest sbottom mass
 * mstau12: squared DR' lightest stau mass
 * mstau22: squared DR' heaviest stau mass
 * msv2   : squared DR' tau sneutrino mass
 *
 * sxt    : sine of DR' stop mixing angle in the MSSM
 * cxt    : cosine of DR' stop mixing angle in the MSSM
 * sxb    : sine of DR' sbottom mixing angle in the MSSM
 * cxb    : cosine of DR' sbottom mixing angle in the MSSM
 * sintau : sine of DR' stau mixing angle in the MSSM
 * costau : cosine of DR' stau mixing angle in the MSSM
 *
 * g3     : DR' strong gauge coupling g3 in the MSSM
 * mu     : DR' mu-parameter in the MSSM (arXiv:0907.4682)
 * tanb   : DR' tan(beta) = vu/vd in the MSSM
 * cotb   : DR' 1/tan(beta) in the MSSM
 * vev2   : squared DR' vev^2 = (vu^2 + vd^2) in the MSSM
 *
 * include_heavy_higgs: factor multiplying the contribution from the heavy Higgs
 */

namespace himalaya {
namespace mssm_twoloophiggs {

/// 2-loop CP-even Higgs contribution O(at*as)
Eigen::Matrix<double, 2, 2> delta_mh2_2loop_at_as(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double g3,
   int include_heavy_higgs);

/// 2-loop CP-even Higgs contribution O(at^2)
Eigen::Matrix<double, 2, 2> delta_mh2_2loop_at_at(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2,
   int include_heavy_higgs);

/// 2-loop CP-even Higgs contribution O(ab*as)
Eigen::Matrix<double, 2, 2> delta_mh2_2loop_ab_as(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double g3,
   int include_heavy_higgs);

/// 2-loop CP-even Higgs contribution O(atau^2)
Eigen::Matrix<double, 2, 2> delta_mh2_2loop_atau_atau(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2,
   int include_heavy_higgs);

/// 2-loop CP-even Higgs contribution O(ab*atau)
Eigen::Matrix<double, 2, 2> delta_mh2_2loop_ab_atau(
   double mtau2, double mb2,
   double mstau12, double mstau22, double msb12, double msb22,
   double sintau, double costau, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2);

} // namespace mssm_twoloophiggs
} // namespace himalaya

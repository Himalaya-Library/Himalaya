// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <iosfwd>
#include <limits>
#include <Eigen/Core>

namespace himalaya {

using V2 = Eigen::Vector2d;   ///< real 2-vector
using RM22 = Eigen::Matrix2d; ///< real 2x2 matrix
using RM33 = Eigen::Matrix3d; ///< real 3x3 matrix
const double NaN = std::numeric_limits<double>::quiet_NaN();

/**
 *         The Himalaya interface struct
 *
 * All input parameters are expected to be provided in the DR'-bar scheme.
 */
struct Parameters {
   // DR'-bar parameters
   double scale{};                  ///< renormalization scale
   double mu{};                          ///< mu parameter, convention of [Phys.Rept. 117 (1985) 75-263]
   double g1{};                          ///< GUT-normalized gauge coupling g1, with gY = g1*Sqrt[3/5]
   double g2{};                          ///< gauge coupling g2
   double g3{};                          ///< gauge coupling g3 SU(3)
   double vd{};                          ///< VEV of down Higgs, with v = Sqrt[vu^2 + vd^2] ~ 246 GeV
   double vu{};                          ///< VEV of up Higgs, with v = Sqrt[vu^2 + vd^2] ~ 246 GeV
   RM33 mq2{RM33::Zero()};          ///< soft-breaking squared left-handed squark mass parameters
   RM33 md2{RM33::Zero()};          ///< soft-breaking squared right-handed down-squark mass parameters
   RM33 mu2{RM33::Zero()};          ///< soft-breaking squared right-handed up-squark mass parameters
   RM33 ml2{RM33::Zero()};          ///< soft-breaking squared left-handed slepton mass parameters
   RM33 me2{RM33::Zero()};          ///< soft-breaking squared right-handed slepton mass parameters
   RM33 Au{RM33::Zero()};          ///< trilinear up type squark-Higgs coupling matrix
   RM33 Ad{RM33::Zero()};          ///< trilinear down type squark-Higgs coupling matrix
   RM33 Ae{RM33::Zero()};          ///< trilinear electron type squark-Higgs coupling matrix
   RM33 Yu{RM33::Zero()};          ///< up-type yukawa coupling matrix
   RM33 Yd{RM33::Zero()};          ///< down-type yukawa coupling matrix
   RM33 Ye{RM33::Zero()};          ///< electron-type yukawa coupling matrix

   // DR'-bar masses
   double M1{};                          ///< bino
   double M2{};                          ///< wino
   double MG{};                          ///< gluino
   double MW{NaN};                ///< W
   double MZ{NaN};                  ///< Z
   double Mt{NaN};                ///< top-quark
   double Mb{NaN};                ///< down-quark
   double Mtau{NaN};                  ///< tau lepton
   double MA{};                          ///< CP-odd Higgs
   V2 MSt{NaN, NaN};              ///< stops
   V2 MSb{NaN, NaN};              ///< sbottoms

   // DR'-bar mixing angles
   double s2t{NaN};               ///< sine of 2 times the stop mixing angle
   double s2b{NaN};               ///< sine of 2 times the sbottom mixing angle

   int massLimit3LThreshold{};           ///< an integer flag to set the mass limit

   double calculateMsq2() const;  ///< calculates average light squark mass squared
   void validate(bool verbose);   ///< validates the parameter set
};

/// prints the Parameters struct to a stream
std::ostream& operator<<(std::ostream&, const Parameters&);

}        //        himalaya

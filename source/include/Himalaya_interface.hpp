#pragma once

#include <iosfwd>
#include <limits>
#include <Eigen/Core>

namespace himalaya {

typedef Eigen::Vector2d V2;   ///< real 2-vector
typedef Eigen::Matrix2d RM22; ///< real 2x2 matrix
typedef Eigen::Matrix3d RM33; ///< real 3x3 matrix

/**
 * 	The Himalaya interface struct
 */
struct Parameters {
   // DR-bar-prime parameters
   double scale{};				/**< renormalization scale */
   double mu{};					/**< mu parameter, convention of [Phys.Rept. 117 (1985) 75-263] */
   double g1{};					/**< GUT-normalized gauge coupling g1, with gY = g1*Sqrt[3/5] */
   double g2{};					/**< gauge coupling g2 */
   double g3{};					/**< gauge coupling g3 SU(3) */
   double vd{};					/**< VEV of down Higgs, with v = Sqrt[vu^2 + vd^2] ~ 246 GeV */
   double vu{};					/**< VEV of up Higgs, with v = Sqrt[vu^2 + vd^2] ~ 246 GeV */
   RM33 mq2{RM33::Zero()};			/**< soft-breaking squared left-handed squark mass parameters */
   RM33 md2{RM33::Zero()};			/**< soft-breaking squared right-handed down-squark mass parameters */
   RM33 mu2{RM33::Zero()};			/**< soft-breaking squared right-handed up-squark mass parameters */
   RM33 ml2{RM33::Zero()};			/**< soft-breaking squared left-handed slepton mass parameters */
   RM33 me2{RM33::Zero()};			/**< soft-breaking squared right-handed slepton mass parameters */
   RM33 Au{RM33::Zero()};			/**< trilinear up type squark-Higgs coupling matrix */
   RM33 Ad{RM33::Zero()};			/**< trilinear down type squark-Higgs coupling matrix */
   RM33 Ae{RM33::Zero()};			/**< trilinear electron type squark-Higgs coupling matrix */
   RM33 Yu{RM33::Zero()};			/**< up-type yukawa coupling matrix */
   RM33 Yd{RM33::Zero()};			/**< down-type yukawa coupling matrix */
   RM33 Ye{RM33::Zero()};			/**< electron-type yukawa coupling matrix */

   // DR-bar-prime masses
   double M1{};					/**< bino */
   double M2{};					/**< wino */
   double MG{};					/**< gluino */
   double MW{std::numeric_limits<double>::quiet_NaN()};		/**< W */
   double MZ{std::numeric_limits<double>::quiet_NaN()};		/**< Z */
   double Mh{std::numeric_limits<double>::quiet_NaN()};		/**< light CP-even Higgs */
   double Mt{std::numeric_limits<double>::quiet_NaN()};		/**< top-quark */
   double Mb{std::numeric_limits<double>::quiet_NaN()};		/**< down-quark */
   double Mtau{std::numeric_limits<double>::quiet_NaN()};	/**< tau lepton */
   double MA{};					/**< CP-odd Higgs */
   V2 MSt{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};	/**< stops */
   V2 MSb{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};	/**< sbottoms */

   // DR-bar-prime mixing angles
   double s2t = std::numeric_limits<double>::quiet_NaN();		/**< sine of 2 times the stop mixing angle */
   double s2b = std::numeric_limits<double>::quiet_NaN();		/**< sine of 2 times the sbottom mixing angle */

   // integer flag to set mass limit
   int massLimit3LThreshold{}; 			/**< an integer flag to set the mass limit */

   double calculateMsq2() const;                /**< calculates average light squark mass squared */
   void validate(bool verbose);                 /**< validates the parameter set */
};

/// prints the Parameters struct to a stream
std::ostream& operator<<(std::ostream&, const Parameters&);

}	//	himalaya

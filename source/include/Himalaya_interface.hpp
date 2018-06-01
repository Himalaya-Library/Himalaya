#pragma once

#include <iosfwd>
#include <limits>
#include <Eigen/Core>

namespace himalaya {

typedef Eigen::Matrix<double,2,1> V2;   ///< real 2-vector
typedef Eigen::Matrix<double,2,2> RM22; ///< real 2x2 matrix
typedef Eigen::Matrix<double,3,3> RM33; ///< real 3x3 matrix

/**
 * 	The Himalaya interface struct
 */
struct Parameters {
   // DR-bar parameters
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
   double At{};					/**< trilinear stop-Higgs coupling */
   double Ab{};					/**< trilinear sbottom-Higgs coupling */
   double Atau{};				/**< trilinear stau-Higgs coupling */

   // DR-bar masses
   double M1{};					/**< bino */
   double M2{};					/**< wino */
   double MG{};					/**< gluino */
   double MW{};					/**< W */
   double MZ{};					/**< Z */
   double Mt{};					/**< top-quark */
   double Mb{};					/**< down-quark */
   double Mtau{};				/**< tau lepton */
   double MA{};					/**< CP-odd Higgs */
   V2 MSt{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};	/**< stops */
   V2 MSb{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};	/**< sbottoms */

   // DR-bar mixing angles
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

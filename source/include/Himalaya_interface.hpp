#pragma once

#include <iosfwd>
#include <limits>
#include <Eigen/Core>

namespace himalaya {

typedef Eigen::Matrix<double,2,1> V2;
typedef Eigen::Matrix<double,2,2> RM22;
typedef Eigen::Matrix<double,3,3> RM33;

/**
 * 	The Himalaya interface struct
 */
struct Parameters {
   // DR-bar parameters
   double scale{};				/**< renormalization scale */
   double mu{};					/**< mu parameter */
   double g3{};					/**< gauge coupling g3 SU(3) */
   double vd{};					/**< VEV of down Higgs */
   double vu{};					/**< VEV of up Higgs */
   RM33 mq2{RM33::Zero()};			/**< soft-breaking squared left-handed squark mass parameters */
   RM33 md2{RM33::Zero()};			/**< soft-breaking squared right-handed down-squark mass parameters */
   RM33 mu2{RM33::Zero()};			/**< soft-breaking squared right-handed up-squark mass parameters */
   double At{};					/**< trilinear stop-Higgs coupling */
   double Ab{};					/**< trilinear sbottom-Higgs coupling */

   // DR-bar masses
   double MG{};					/**< gluino */
   double MW{};					/**< W */
   double MZ{};					/**< Z */
   double Mt{};					/**< top-quark */
   double Mb{};					/**< down-quark */
   double MA{};					/**< CP-odd Higgs */
   V2 MSt{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};	/**< stops */
   V2 MSb{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};	/**< sbottoms */

   // DR-bar mixing angles
   double s2t = std::numeric_limits<double>::quiet_NaN();		/**< sine of 2 times the stop mixing angle */
   double s2b = std::numeric_limits<double>::quiet_NaN();		/**< sine of 2 times the sbottom mixing angle */

   // integer flag to set mass limit
   int massLimit3LThreshold{}; 			/**< an integer flag to set the mass limit */

   double calculateMsq2() const;                /** calculates average light squark mass squared */
   void validate(bool verbose);                 /** validates the parameter set */
};

std::ostream& operator<<(std::ostream&, const Parameters&);

}	//	himalaya

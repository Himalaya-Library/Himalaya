#ifndef H3m_interface_HPP
#define H3m_interface_HPP

#include <complex>
#include <Eigen>

namespace h3m {

typedef Eigen::Matrix<double,2,1> V2;
typedef Eigen::Matrix<double,2,2> RM22;
typedef Eigen::Matrix<double,3,3> RM33;

struct Parameters {
   // DR-bar parameters
   double scale{};         ///< renormalization scale
   double mu{};            ///< mu parameter
   double g3{};            ///< gauge coupling g3 SU(3)
   double vd{};            ///< VEV of down Higgs
   double vu{};            ///< VEV of up Higgs
   RM33 mq2{RM33::Zero()}; ///< soft-breaking squared left-handed squark mass parameters
   RM33 md2{RM33::Zero()}; ///< soft-breaking squared right-handed down-squark mass parameters
   RM33 mu2{RM33::Zero()}; ///< soft-breaking squared right-handed up-squark mass parameters
   double At{};
   double Ab{};

   // DR-bar masses
   double MG{};            ///< gluino
   double MW{};            ///< W
   double MZ{};            ///< Z
   double Mt{};            ///< top-quark
   double Mb{};            ///< down-quark
   double MA{};
   V2 MSt{V2::Zero()};     ///< stops
   V2 MSb{V2::Zero()};     ///< sbottoms

   // DR-bar mixing matrices
   RM22 Zt{RM22::Zero()};  ///< stop mixing matrix
   RM22 Zb{RM22::Zero()};  ///< sbottom mixing matrix
};

/// returns 3-loop corrections to CP-even Higgs pole masses in the MSSM
RM22 calculate_DMh_3L(const Parameters& p);

}	//	h3m

#endif	//	H3m_interface_HPP

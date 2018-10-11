// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "Himalaya_interface.hpp"
#include <Eigen/Core>

namespace himalaya {
namespace mh1l {

class MSSM_mass_eigenstates {
public:
   MSSM_mass_eigenstates(const Parameters&);

   /// Higgs 1-loop contribution DR'
   RM22 delta_mh2_1loop(double p2) const;
   /// Higgs 1-loop contribution DR' for p = g1 = g2 = 0
   RM22 delta_mh2_1loop_effpot_gaugeless() const;

private:
   using A2 = Eigen::Array<double,2,1>;

   Parameters pars;         ///< MSSM DR' parameters

   double MFt{0.};          ///< MSSM DR' top mass
   double MFb{0.};          ///< MSSM DR' bottom mass
   double MFtau{0.};        ///< MSSM DR' tau mass
   A2 M2St{V2::Zero()};     ///< MSSM DR' squared stop masses
   A2 M2Sb{V2::Zero()};     ///< MSSM DR' squared sbottom masses
   A2 M2Stau{V2::Zero()};   ///< MSSM DR' squared stau masses

   RM22 ZT{RM22::Zero()};   ///< MSSM DR' stop mixing matrix
   RM22 ZB{RM22::Zero()};   ///< MSSM DR' sbottom mixing matrix
   RM22 ZTau{RM22::Zero()}; ///< MSSM DR' stau mixing matrix

   /// calculates all DR' masses and mixings
   void calculate_parameters();

   void calculate_MFt();    ///< calculates DR' top mass
   void calculate_MFb();    ///< calculates DR' bottom mass
   void calculate_MFtau();  ///< calculates DR' tau mass
   void calculate_MSt();    ///< calculates DR' stop mass
   void calculate_MSb();    ///< calculates DR' sbottom mass
   void calculate_MStau();  ///< calculates DR' stau mass

   RM22 get_mass_matrix_St() const;   ///< stop mass matrix
   RM22 get_mass_matrix_Sb() const;   ///< sbottom mass matrix
   RM22 get_mass_matrix_Stau() const; ///< stau mass matrix

   /// A0 Passarino-Veltman function
   double A0(double) const;
   /// B0 Passarino-Veltman function
   double B0(double, double, double) const;
};

} // namespace mh1l
} // namespace himalaya

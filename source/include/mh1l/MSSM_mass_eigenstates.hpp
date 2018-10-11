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
   using A4 = Eigen::Array<double,4,1>;

   Parameters pars;         ///< MSSM DR' parameters

   double M2VWm{0.};        ///< MSSM DR' squared W mass
   double M2VZ{0.};         ///< MSSM DR' squared Z mass
   double MFt{0.};          ///< MSSM DR' top mass
   double MFb{0.};          ///< MSSM DR' bottom mass
   double MFtau{0.};        ///< MSSM DR' tau mass
   double M2SveL{0.};       ///< MSSM DR' electron-like sneutrino mass
   double M2SvmL{0.};       ///< MSSM DR' muon-like sneutrino mass
   double M2SvtL{0.};       ///< MSSM DR' tau-like sneutrino mass
   A2 M2Su{V2::Zero()};     ///< MSSM DR' squared sup masses
   A2 M2Sd{V2::Zero()};     ///< MSSM DR' squared sdown masses
   A2 M2Sc{V2::Zero()};     ///< MSSM DR' squared scharm masses
   A2 M2Ss{V2::Zero()};     ///< MSSM DR' squared sstrange masses
   A2 M2St{V2::Zero()};     ///< MSSM DR' squared stop masses
   A2 M2Sb{V2::Zero()};     ///< MSSM DR' squared sbottom masses
   A2 M2Se{V2::Zero()};     ///< MSSM DR' squared selectron masses
   A2 M2Sm{V2::Zero()};     ///< MSSM DR' squared smuon masses
   A2 M2Stau{V2::Zero()};   ///< MSSM DR' squared stau masses
   A2 M2hh{A2::Zero()};     ///< MSSM DR' CP-even higgs mass
   A2 M2Ah{A2::Zero()};     ///< MSSM DR' CP-odd higgs mass
   A2 M2Hpm{A2::Zero()};    ///< MSSM DR' charged higgs mass
   A4 MChi{A4::Zero()};     ///< MSSM DR' neutralino mass
   A2 MCha{A2::Zero()};     ///< MSSM DR' chargino mass

   RM22 ZU{RM22::Zero()};   ///< MSSM DR' sup mixing matrix
   RM22 ZD{RM22::Zero()};   ///< MSSM DR' sdown mixing matrix
   RM22 ZC{RM22::Zero()};   ///< MSSM DR' scharm mixing matrix
   RM22 ZS{RM22::Zero()};   ///< MSSM DR' sstrange mixing matrix
   RM22 ZT{RM22::Zero()};   ///< MSSM DR' stop mixing matrix
   RM22 ZB{RM22::Zero()};   ///< MSSM DR' sbottom mixing matrix
   RM22 ZE{RM22::Zero()};   ///< MSSM DR' selectron mixing matrix
   RM22 ZM{RM22::Zero()};   ///< MSSM DR' smuon mixing matrix
   RM22 ZTau{RM22::Zero()}; ///< MSSM DR' stau mixing matrix
   RM22 ZH{RM22::Zero()};   ///< MSSM DR' CP-even Higgs mixing matrix
   RM22 ZA{RM22::Zero()};   ///< MSSM DR' CP-odd Higgs mixing matrix
   RM22 ZP{RM22::Zero()};   ///< MSSM DR' charged Higgs mixing matrix
   CM44 ZN{CM44::Zero()};   ///< MSSM DR' neutralino mixing matrix
   CM22 UM{CM22::Zero()};   ///< MSSM DR' positive chargino mixing matrix
   CM22 UP{CM22::Zero()};   ///< MSSM DR' negative chargino mixing matrix

   /// calculates all DR' masses and mixings
   void calculate_parameters();

   void calculate_MVWm();   ///< calculates DR' W mass
   void calculate_MVZ();    ///< calculates DR' Z mass
   void calculate_MFt();    ///< calculates DR' top mass
   void calculate_MFb();    ///< calculates DR' bottom mass
   void calculate_MFtau();  ///< calculates DR' tau mass
   void calculate_MSveL();  ///< calculates DR' electron-like sneutrino mass
   void calculate_MSvmL();  ///< calculates DR' muon-like sneutrino mass
   void calculate_MSvtL();  ///< calculates DR' tau-like sneutrino mass
   void calculate_MSu();    ///< calculates DR' sup masses
   void calculate_MSd();    ///< calculates DR' sdown masses
   void calculate_MSc();    ///< calculates DR' scharm masses
   void calculate_MSs();    ///< calculates DR' sstrange masses
   void calculate_MSt();    ///< calculates DR' stop masses
   void calculate_MSb();    ///< calculates DR' sbottom masses
   void calculate_MSe();    ///< calculates DR' selectron masses
   void calculate_MSm();    ///< calculates DR' smuon masses
   void calculate_MStau();  ///< calculates DR' stau masses
   void calculate_Mhh();    ///< calculates DR' CP-even Higgs masses
   void calculate_MAh();    ///< calculates DR' CP-odd Higgs masses
   void calculate_MHpm();   ///< calculates DR' charged Higgs masses
   void calculate_MChi();   ///< calculates DR' neutralino masses
   void calculate_MCha();   ///< calculates DR' chargino masses

   RM22 get_mass_matrix_Su() const;   ///< sup mass matrix
   RM22 get_mass_matrix_Sd() const;   ///< sdown mass matrix
   RM22 get_mass_matrix_Sc() const;   ///< scharm mass matrix
   RM22 get_mass_matrix_Ss() const;   ///< sstrange mass matrix
   RM22 get_mass_matrix_St() const;   ///< stop mass matrix
   RM22 get_mass_matrix_Sb() const;   ///< sbottom mass matrix
   RM22 get_mass_matrix_Se() const;   ///< selectron mass matrix
   RM22 get_mass_matrix_Sm() const;   ///< smuon mass matrix
   RM22 get_mass_matrix_Stau() const; ///< stau mass matrix
   RM22 get_mass_matrix_hh() const;   ///< CP-even Higgs mass matrix
   RM22 get_mass_matrix_Ah() const;   ///< CP-odd Higgs mass matrix
   RM22 get_mass_matrix_Hpm() const;  ///< charged Higgs mass matrix
   RM44 get_mass_matrix_Chi() const;  ///< neutralino mass matrix
   RM22 get_mass_matrix_Cha() const;  ///< chargino mass matrix

   /// A0 Passarino-Veltman function
   double A0(double) const;
   /// B0 Passarino-Veltman function
   double B0(double, double, double) const;
};

} // namespace mh1l
} // namespace himalaya

// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "Himalaya_interface.hpp"
#include <map>
#include <tuple>
#include <Eigen/Core>

namespace himalaya {
namespace mh1l {

struct MSSM_spectrum {
   using A2 = Eigen::Array<double,2,1>;
   using A4 = Eigen::Array<double,4,1>;

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
   RM44 ZN{RM44::Zero()};   ///< MSSM DR' neutralino mixing matrix
   RM22 UM{RM22::Zero()};   ///< MSSM DR' positive chargino mixing matrix
   RM22 UP{RM22::Zero()};   ///< MSSM DR' negative chargino mixing matrix

   void calculate_spectrum(const Parameters&);

   void calculate_MVWm(const Parameters&);   ///< calculates DR' W mass
   void calculate_MVZ(const Parameters&);    ///< calculates DR' Z mass
   void calculate_MFt(const Parameters&);    ///< calculates DR' top mass
   void calculate_MFb(const Parameters&);    ///< calculates DR' bottom mass
   void calculate_MFtau(const Parameters&);  ///< calculates DR' tau mass
   void calculate_MSveL(const Parameters&);  ///< calculates DR' electron-like sneutrino mass
   void calculate_MSvmL(const Parameters&);  ///< calculates DR' muon-like sneutrino mass
   void calculate_MSvtL(const Parameters&);  ///< calculates DR' tau-like sneutrino mass
   void calculate_MSu(const Parameters&);    ///< calculates DR' sup masses
   void calculate_MSd(const Parameters&);    ///< calculates DR' sdown masses
   void calculate_MSc(const Parameters&);    ///< calculates DR' scharm masses
   void calculate_MSs(const Parameters&);    ///< calculates DR' sstrange masses
   void calculate_MSt(const Parameters&);    ///< calculates DR' stop masses
   void calculate_MSb(const Parameters&);    ///< calculates DR' sbottom masses
   void calculate_MSe(const Parameters&);    ///< calculates DR' selectron masses
   void calculate_MSm(const Parameters&);    ///< calculates DR' smuon masses
   void calculate_MStau(const Parameters&);  ///< calculates DR' stau masses
   void calculate_Mhh(const Parameters&);    ///< calculates DR' CP-even Higgs masses
   void calculate_MAh(const Parameters&);    ///< calculates DR' CP-odd Higgs masses
   void calculate_MHpm(const Parameters&);   ///< calculates DR' charged Higgs masses
   void calculate_MChi(const Parameters&);   ///< calculates DR' neutralino masses
   void calculate_MCha(const Parameters&);   ///< calculates DR' chargino masses

   RM22 get_mass_matrix_Su(const Parameters&) const;   ///< sup mass matrix
   RM22 get_mass_matrix_Sd(const Parameters&) const;   ///< sdown mass matrix
   RM22 get_mass_matrix_Sc(const Parameters&) const;   ///< scharm mass matrix
   RM22 get_mass_matrix_Ss(const Parameters&) const;   ///< sstrange mass matrix
   RM22 get_mass_matrix_St(const Parameters&) const;   ///< stop mass matrix
   RM22 get_mass_matrix_Sb(const Parameters&) const;   ///< sbottom mass matrix
   RM22 get_mass_matrix_Se(const Parameters&) const;   ///< selectron mass matrix
   RM22 get_mass_matrix_Sm(const Parameters&) const;   ///< smuon mass matrix
   RM22 get_mass_matrix_Stau(const Parameters&) const; ///< stau mass matrix
   RM22 get_mass_matrix_hh(const Parameters&) const;   ///< CP-even Higgs mass matrix
   RM22 get_mass_matrix_Ah(const Parameters&) const;   ///< CP-odd Higgs mass matrix
   RM22 get_mass_matrix_Hpm(const Parameters&) const;  ///< charged Higgs mass matrix
   RM44 get_mass_matrix_Chi(const Parameters&) const;  ///< neutralino mass matrix
   RM22 get_mass_matrix_Cha(const Parameters&) const;  ///< chargino mass matrix
};

/// prints the spectrum
std::ostream& operator<<(std::ostream&, const MSSM_spectrum&);

class MSSM_mass_eigenstates {
public:
   MSSM_mass_eigenstates(const Parameters&);

   /// calculates squared Higgs masses
   std::tuple<V2,V2,V2> calculate_Mh2() const;
   /// Higgs 1-loop contribution DR'
   RM22 delta_mh2_1loop(double p2) const;
   /// Higgs 1-loop contribution DR' for p = g1 = g2 = 0
   RM22 delta_mh2_1loop_gaugeless() const;
   /// derivative of Higgs 1-loop contribution DR' for p = g1 = g2 = 0
   RM22 delta_mh2_1loop_gaugeless_deriv() const;
   /// Higgs 2-loop contributions DR' for p = 0
   RM22 delta_mh2_2loop() const;
   /// Higgs 2-loop contributions DR' for p = g1 = g2 = 0 from momentum iteration
   RM22 delta_mh2_2loop_mom_it() const;
   /// returns CP-even Higgs mass matrix
   RM22 get_mass_matrix_hh() const;
   /// enable/disable loop corrections
   void set_correction(int, int);
   /// enable/disable momentum iteration
   void enable_mom_it(bool);

   friend std::ostream& operator<<(std::ostream&, const MSSM_mass_eigenstates&);

private:
   Parameters pars;            ///< MSSM DR' parameters
   MSSM_spectrum masses;       ///< MSSM DR' masses / mixings
   MSSM_spectrum gaugeless;    ///< MSSM DR' masses / mixings for g1 = g2 = 0
   std::map<int,int> orders{}; ///< enable/disable corrections
   bool include_mom_it{true};  ///< include/exclude momentum iteration

   /// calculates tree-level squared Higgs masses
   V2 calculate_Mh2_tree() const;

   /// calculates all DR' masses and mixings
   void calculate_parameters();

   /// sets g1 = g2 = 0
   Parameters make_gaugeless(const Parameters&) const;

   /// A0 Passarino-Veltman function
   double A0(double) const;
   /// B0 Passarino-Veltman function
   double B0(double, double, double) const;
   /// derivative of B0 function w.r.t. p^2, for p^2 = 0
   double D1B0(double, double) const;
};

/// prints the internals of MSSM_mass_eigenstates
std::ostream& operator<<(std::ostream&, const MSSM_mass_eigenstates&);

} // namespace mh1l
} // namespace himalaya

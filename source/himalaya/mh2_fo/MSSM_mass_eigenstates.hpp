// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "himalaya/Himalaya_interface.hpp"

#include "himalaya/mh2_fo/MSSM_spectrum.hpp"
#include "himalaya/misc/CouplingOrders.hpp"

#include <array>
#include <iosfwd>
#include <tuple>
#include <Eigen/Core>

/**
 * @file MSSM_mass_eigenstates.hpp
 *
 * @brief Contains the definition of the \a
 * MSSM_mass_eigenstates class.
 */

namespace himalaya {
namespace mh2_fo {

/// momentum iteration settings
enum class Momentum_iteration {
   off,  ///< no momentum iteration
   pert, ///< perturbatively up to 2-loop level
   num,  ///< numerically (all orders)
};

/// diagonalization settings
enum class Diagonalization {
   pert, ///< perturbatively
   num,  ///< numerically
};

/**
 * @class MSSM_mass_eigenstates
 *
 * @brief This class performs a fixed-order calculation of the light
 * CP-even Higgs mass up to 2-loop order.
 */
class MSSM_mass_eigenstates {
public:
   MSSM_mass_eigenstates(const Parameters&, bool only_at = false);

   /// calculates squared Higgs masses
   std::tuple<double,double,double> calculate_Mh2() const;
   /// returns tree-level CP-even Higgs mass matrix
   RM22 get_mass_matrix_hh() const;
   /// returns tree-level CP-even Higgs mass matrix for p = g1 = g2 = 0
   RM22 get_mass_matrix_hh_gaugeless() const;
   /// Higgs 1-loop contribution DR'
   RM22 delta_mh2_1loop(double p2) const;
   /// Higgs 1-loop contribution DR' for p = g1 = g2 = 0
   RM22 delta_mh2_1loop_gaugeless() const;
   /// derivative of Higgs 1-loop contribution DR' for p = g1 = g2 = 0
   RM22 delta_mh2_1loop_gaugeless_deriv() const;
   /// Higgs 2-loop contributions DR' for p = 0
   RM22 delta_mh2_2loop() const;
   /// Higgs 2-loop contributions DR' for p = g1 = g2 = 0 from momentum iteration
   RM22 delta_mh2_2loop_mom_it_pert() const;
   /// Higgs 2-loop (and higher) contributions DR' from numerical momentum iteration
   RM22 delta_mh2_2loop_mom_it_num(double precision_goal = 1e-5, int max_iterations = 100) const;
   /// enable/disable loop corrections
   void set_correction(CouplingOrders::CouplingOrders, int);
   /// customize diagonalization
   void set_diagonalization(Diagonalization);
   /// customize momentum iteration
   void set_mom_it(Momentum_iteration, double mom_it_precision_goal_ = 1e-5, int mom_it_max_iterations_ = 100);

   friend std::ostream& operator<<(std::ostream&, const MSSM_mass_eigenstates&);

private:
   Parameters pars;            ///< MSSM DR' parameters
   MSSM_spectrum masses;       ///< MSSM DR' masses / mixings
   MSSM_spectrum gaugeless;    ///< MSSM DR' masses / mixings for g1 = g2 = 0
   std::array<int,CouplingOrders::NUMBER_OF_COUPLING_ORDERS> orders{};   ///< enable/disable corrections
   Momentum_iteration mom_it{Momentum_iteration::pert}; ///< momentum iteration settings
   double mom_it_precision_goal{1e-5}; ///< precision goal for numeric momentum iteration
   int mom_it_max_iterations{100};     ///< maximum number of numeric momentum iterations
   Diagonalization diagonalization{Diagonalization::pert}; ///< diagonalization settings

   /// calculates tree-level squared Higgs masses
   V2 calculate_Mh2_tree() const;

   /// sets g1 = g2 = 0
   static Parameters make_gaugeless(const Parameters&);

   /// sets yb = ytau = 0
   static Parameters make_3rd_gen(const Parameters&);

   /// A0 Passarino-Veltman function
   double A0(double) const;
   /// B0 Passarino-Veltman function
   double B0(double, double, double) const;
   /// derivative of B0 function w.r.t. p^2, for p^2 = 0
   double D1B0(double, double) const;
};

/// prints the internals of MSSM_mass_eigenstates
std::ostream& operator<<(std::ostream&, const MSSM_mass_eigenstates&);

} // namespace mh2_fo
} // namespace himalaya

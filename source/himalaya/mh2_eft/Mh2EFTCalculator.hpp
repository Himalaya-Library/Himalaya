// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "himalaya/Himalaya_interface.hpp"

#include "himalaya/misc/CouplingOrders.hpp"

#include <iosfwd>
#include <array>

/**
 * @file Mh2EFTCalculator.hpp
 * @brief Definition of EFT Higgs mass calculation class.
 */

namespace himalaya{
namespace mh2_eft{

/**
 * The Mh2 EFT calculator class
 */
class Mh2EFTCalculator{
public:
   /**
    * Constructor
    * @param p_ a HimalayaInterface struct
    * @param verbose a bool enable the output of the parameter validation. Enabled by default
    */
   Mh2EFTCalculator(const Parameters& p_, bool verbose = true);
   /**
    * Returns the tree-level EFT contribution to the Higgs mass
    */
   double getDeltaMh2EFT0Loop() const;
   /**
    * Returns the 1-loop EFT contribution to the Higgs mass
    * @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
    * @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
    */
   double getDeltaMh2EFT1Loop(int omitSMLogs, int omitMSSMLogs) const;
   /**
    * Returns the 1-loop EFT contribution to the Higgs mass
    * @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
    * @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
    */
   double getDeltaMh2EFT2Loop(int omitSMLogs, int omitMSSMLogs) const;
   /**
    * Returns the 1-loop EFT contribution to the Higgs mass
    * @param omitSMLogs an integer flag to remove all Log(mu^2/mt^2) terms
    * @param omitMSSMLogs an integer flag to remove all Log(mu^2/Mx^2) terms
    * @param omitDeltaLambda3L an integer flag to disable the MSSM contribution to delta_lambda_3L
    */
   double getDeltaMh2EFT3Loop(int omitSMLogs, int omitMSSMLogs, int omitDeltaLambda3L = 1) const;
   /**
    * Returns the matching relation of delta_lambda 3L for the degenerate mass case
    * @param scale the renormalization scale
    * @param mst1 the mass of the light stop quark
    * @param Xt stop mixing parameter
    * @param omitlogs an integer flag to remove all logarithmic terms
    * @return delta_Lambda 3L
    */
   double getDeltaLambdaDegenerate(double scale, double mst1, double Xt, int omitlogs) const;
   /**
    * Sets the flag to enable or disable a correction of a given variable
    * @param order an integer taken from the CouplingOrders enum
    * @param flag set to 1 to enable and to 0 to disable the chosen correction
    */
   void setCorrectionFlag(CouplingOrders::CouplingOrders order, int flag);

   friend std::ostream& operator<<(std::ostream&, const Mh2EFTCalculator&);
private:
   /// calculate beta
   double calcBeta() const;
   /// calculate tan(beta)
   double calcTanBeta() const;
   /// calculate v
   double calcV() const;
   /// calculate v^2
   double calcV2() const;
   /// calculate sin(beta)
   double calcSinBeta() const;
   /// calculate cos(beta)
   double calcCosBeta() const;

   Parameters p{}; ///< The HimalayaInterface struct
   std::array<int, CouplingOrders::NUMBER_OF_COUPLING_ORDERS> orders{}; ///< holds all CouplingOrders to enable/disable certain corrections
};

/// prints loop corrections for v^2 << MS^2
std::ostream& operator<<(std::ostream&, const Mh2EFTCalculator&);

} // namespace mh2_eft
} // namespace himalaya

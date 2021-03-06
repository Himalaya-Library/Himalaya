// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "himalaya/Himalaya_interface.hpp"
#include "himalaya/HierarchyObject.hpp"
#include <array>
#include <Eigen/Core>

/**
 * @file HierarchyCalculator.hpp
 *
 * @brief Definition of the HierarchyCalculator.
 *
 * Definition of the HierarchyCalculatur, which selects the best
 * suited hierarchy for the given MSSM parameter point and calculates
 * the loop corrections.
 */

namespace himalaya {
   /**
    * The HierarchyCalculatur class
    */
   class HierarchyCalculator{
   public:
      /**
       * Constructor
       * @param p_ Himalaya input parameters
       * @param verbose_ suppress informative output during the calculation, if set to false
       */
      HierarchyCalculator(const Parameters& p_, bool verbose_ = true);
      /**
       * Calculates the 3-loop mass matrix and other information of the hierarchy selection process.
       * @param isAlphab a bool which determines if the returned object is proportinal to alpha_b.
       * @return A HierarchyObject which holds all information of the calculation.
       */
      HierarchyObject calculateDMh3L(bool isAlphab);
      /**
       * Compares deviation of all hierarchies with the exact two-loop result and returns the hierarchy which minimizes the error.
       * @param ho a HierarchyObject with constant isAlphab.
       * @return An integer which is identified with the suitable hierarchy.
       */
      int compareHierarchies(HierarchyObject& ho);
      /**
       * Calculates the hierarchy contributions for a specific hierarchy at a specific loop order.
       * @param ho a HierarchyObject with constant isAlphab.
       * @param oneLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded one-loop results to the returned value, respectivley.
       * @param twoLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded two-loop results to the returned value, respectivley.
       * @param threeLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded three-loop results to the returned value, respectivley.
       * @throws runtime_error Throws a runtime_error if the tree-level is requested in terms of hierarchies.
       * @return The loop corrected Higgs mass matrix which contains the expanded corrections at the given order.
       */
      Eigen::Matrix2d calculateHierarchy(himalaya::HierarchyObject& ho, int oneLoopFlagIn, int twoLoopFlagIn, int threeLoopFlagIn) const;
      /**
       * Calculates the contribution to the order (alpha_x) and (alpha_s alpha_x) as the difference
       * of the Higgs mass matrices of the MDR and DR scheme. Here, x can be t or b.
       * @param ho a HierarchyObject with constant isAlphab.
       * @param shiftOneLoop a bool to shift the terms at one-loop level.
       * @param shiftTwoLoop a bool to shift the terms at two-loop level.
       * @return The loop corrected Higgs mass matrix difference of the MDR and DR scheme at the given order.
       */
      Eigen::Matrix2d calcDRbarToMDRbarShift(const HierarchyObject& ho, bool shiftOneLoop, bool shiftTwoLoop) const;
      /**
       * Calculates the loop corrected Higgs mass matrix at the order O(alpha_x). Here, x can be t or b.
       * @param ho a HierarchyObject with constant isAlphab.
       * @param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
       * @param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
       * @return The loop corrected Higgs mass matrix at the order O(alpha_x).
       */
      Eigen::Matrix2d getMt41L(const HierarchyObject& ho, unsigned shiftOneLoop, unsigned shiftTwoLoop) const;
      /**
       * Calculates the loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s). Here, x can be t or b.
       * @param ho a HierarchyObject with constant isAlphab.
       * @param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
       * @param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
       * @return The loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s).
       */
      Eigen::Matrix2d getMt42L(const HierarchyObject& ho, unsigned shiftOneLoop, unsigned shiftTwoLoop) const;
      /**
       * Shifts Msx1 according to the hierarchy to the MDR scheme.
       * @param ho a HierarchyObject with constant isAlphab.
       * @param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
       * @param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
       * @return A double which is the MDR sx_1 mass.
       */
      double shiftMst1ToMDR(const HierarchyObject& ho, unsigned oneLoopFlag, unsigned twoLoopFlag) const;
      /**
       * Shifts Mst2 according to the hierarchy to the MDR scheme.
       * @param ho a HierarchyObject with constant isAlphab.
       * @param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
       * @param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
       * @return A double which is the MDR stop_2 mass.
       */
      double shiftMst2ToMDR(const HierarchyObject& ho, unsigned oneLoopFlag, unsigned twoLoopFlag) const;
      /**
       * Estimates the uncertainty of the expansion at a given order.
       * @param ho a HierarchyObject with constant isAlphab.
       * @param massMatrix the CP-even Higgs mass matrix without the corrections whose uncertainty should be estimated.
       * @param oneLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the one-loop expansion terms.
       * @param twoLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the two-loop expansion terms.
       * @param threeLoopFlag an integer flag which is 0 or 1 in order to estimte the uncertainty of the three-loop expansion terms.
       * @return A double which is the estimated uncertainty.
       */
      double getExpansionUncertainty(const himalaya::HierarchyObject& ho, const Eigen::Matrix2d& massMatrix,
                                     int oneLoopFlag, int twoLoopFlag, int threeLoopFlag);
   private:
      Parameters p{};                    ///< Himalaya input parameters
      bool verbose{true};                ///< enable/disable verbose output
      std::array<int, 12> expansionDepth{}; ///< hierarchy expansion depth
      /**
       * Checks if a hierarchy is suitable to the given mass spectrum.
       * @param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
       * @returns A bool if the hierarchy candidate is suitable to the given mass spectrum.
       */
      bool isHierarchySuitable(const HierarchyObject& ho) const;
      /**
       * Shifts the H3m renormalization scheme to DR' scheme
       * @param ho a HierarchyObject with constant isAlphab
       * @return A matrix which shifts the H3m scheme to the DR' scheme at three-loop level
       *
       */
      Eigen::Matrix2d shiftH3mToDRbarPrime(const HierarchyObject& ho) const;
      /**
       * Shifts the H3m renormalization scheme to DR' scheme
       * @param ho a HierarchyObject with constant isAlphab
       * @param omitLogs a flag to omit logarithmic contributions: (0) omit logs, (1) add them
       * @return A double which shifts the H3m scheme to the DR' scheme at three-loop level
       *
       */
      double shiftH3mToDRbarPrimeMh2(const himalaya::HierarchyObject& ho, int omitLogs) const;
      /**
       * Fills in delta_lambda @ 3L to the given HierarchyObject
       * @param ho a HierrachyObject
       * @param omitXtOrders a bool to omit xtOrders of delta_lambda_EFT
       */
      void calcDeltaLambda3L(HierarchyObject& ho, bool omitXtOrders) const;
      /// calculate sin(beta)
      double calcSinBeta() const;
      /// calculate cos(beta)
      double calcCosBeta() const;
      /// calculate tan(beta)
      double calcTanBeta() const;
      /// calculate beta from tan(beta)
      double calcBeta() const;
      /// calculate v = sqrt(vu^2 + vd^2)
      double calcV() const;
      /// calculate v^2 = vu^2 + vd^2
      double calcV2() const;
      /// calculate prefactor of the Higgs mass matrix
      double calcHiggsMassMatrixPrefactor() const;
      /// calculate prefactor as/(4 Pi)
      double calcAsOver4Pi() const;
      /// mean (non-squared) light squark mass
      double calcMeanMsq() const;
      /// calculate sfermion masses shifted to MDR
      std::array<double, 2> calcMsfMDRFlag(const HierarchyObject& ho, int loopOrder) const;
  };
}        // himalaya

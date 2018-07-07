// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "Himalaya_interface.hpp"
#include "HierarchyObject.hpp"
#include "version.hpp"
#include <map>

namespace himalaya{
   /**
    * The HierarchyCalculatur class
    */
   class HierarchyCalculator{
   public:
      /**
       *         Constructor
       *         @param p_ Himalaya input parameters
       *         @param verbose_ suppress informative output during the calculation, if set to false
       */
      HierarchyCalculator(const Parameters& p_, const bool verbose_ = true);
      /**
       *         Calculates the 3-loop mass matrix and other information of the hierarchy selection process.
       *         @param isAlphab a bool which determines if the returned object is proportinal to alpha_b.
       *         @return A HierarchyObject which holds all information of the calculation.
       */
      HierarchyObject calculateDMh3L(bool isAlphab);
      /**
       *         Compares deviation of all hierarchies with the exact two-loop result and returns the hierarchy which minimizes the error.
       *         @param ho a HierarchyObject with constant isAlphab.
       *         @return An integer which is identified with the suitable hierarchy.
       */
      int compareHierarchies(HierarchyObject& ho);
      /**
       *         Calculates the hierarchy contributions for a specific hierarchy at a specific loop order.
       *         @param ho a HierarchyObject with constant isAlphab.
       *         @param oneLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded one-loop results to the returned value, respectivley.
       *         @param twoLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded two-loop results to the returned value, respectivley.
       *         @param threeLoopFlagIn an integer flag which is 0 or 1 in order to add or omit the expanded three-loop results to the returned value, respectivley.
       *         @throws runtime_error Throws a runtime_error if the tree-level is requested in terms of hierarchies.
       *         @return The loop corrected Higgs mass matrix which contains the expanded corrections at the given order.
       */
      Eigen::Matrix2d calculateHierarchy(himalaya::HierarchyObject& ho, const int oneLoopFlagIn, const int twoLoopFlagIn, const int threeLoopFlagIn) const;
      /**
       *         Calculates the contribution to the order (alpha_x) and (alpha_s alpha_x) as the difference
       *        of the Higgs mass matrices of the MDR and DR scheme. Here, x can be t or b.
       *         @param ho a HierarchyObject with constant isAlphab.
       *         @param shiftOneLoop a bool to shift the terms at one-loop level.
       *         @param shiftTwoLoop a bool to shift the terms at two-loop level.
       *         @return The loop corrected Higgs mass matrix difference of the MDR and DR scheme at the given order.
       */
      Eigen::Matrix2d calcDRbarToMDRbarShift(const HierarchyObject& ho, const bool shiftOneLoop, const bool shiftTwoLoop) const;
      /**
       *         Calculates the loop corrected Higgs mass matrix at the order O(alpha_x). Here, x can be t or b.
       *         @param ho a HierarchyObject with constant isAlphab.
       *         @param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
       *         @param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
       *         @return The loop corrected Higgs mass matrix at the order O(alpha_x).
       */
      Eigen::Matrix2d getMt41L(const HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop) const;
      /**
       *         Calculates the loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s). Here, x can be t or b.
       *         @param ho a HierarchyObject with constant isAlphab.
       *         @param shiftOneLoop An integer flag which is 0 or 1 in order to shift the one-loop terms to the MDR scheme.
       *         @param shiftTwoLoop An integer flag which is 0 or 1 in order to shift the two-loop terms to the MDR scheme.
       *         @return The loop corrected Higgs mass matrix at the order O(alpha_x*alpha_s).
       */
      Eigen::Matrix2d getMt42L(const HierarchyObject& ho, const unsigned int shiftOneLoop, const unsigned int shiftTwoLoop) const;
      /**
       *        Shifts Msx1 according to the hierarchy to the MDR scheme.
       *         @param ho a HierarchyObject with constant isAlphab.
       *         @param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
       *         @param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
       *         @return A double which is the MDR sx_1 mass.
       */
      double shiftMst1ToMDR(const HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag) const;
      /**
       *        Shifts Mst2 according to the hierarchy to the MDR scheme.
       *         @param ho a HierarchyObject with constant isAlphab.
       *         @param oneLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s).
       *         @param twoLoopFlag an integer flag which is 0 or 1 in order to shift the order O(alpha_s^2).
       *         @return A double which is the MDR stop_2 mass.
       */
      double shiftMst2ToMDR(const HierarchyObject& ho, const unsigned int oneLoopFlag, const unsigned int twoLoopFlag) const;
      /**
       *         Estimates the uncertainty of the expansion at a given order.
       *         @param ho a HierarchyObject with constant isAlphab.
       *         @param massMatrix the CP-even Higgs mass matrix without the corrections whose uncertainty should be estimated.
       *         @param oneLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the one-loop expansion terms.
       *         @param twoLoopFlag an integer flag which is 0 or 1 in order to estimate the uncertainty of the two-loop expansion terms.
       *         @param threeLoopFlag an integer flag which is 0 or 1 in order to estimte the uncertainty of the three-loop expansion terms.
       *         @return A double which is the estimated uncertainty.
       */
      double getExpansionUncertainty(himalaya::HierarchyObject& ho, const Eigen::Matrix2d& massMatrix, const unsigned int oneLoopFlag,
                                     const unsigned int twoLoopFlag, const unsigned int threeLoopFlag);
   private:
      Parameters p{};     ///< Himalaya input parameters
      double Al4p{};      ///< alpha_s/(4*Pi)
      double lmMgl{};     ///< log(pow2(p.scale / Mgl))
      double lmMsq{};     ///< log(pow2(p.scale / Msq))
      double Mgl{};       ///< Gluino mass
      double Msq{};       ///< mean light squark mass
      double prefac{};    ///< prefactor of the Higgs mass matrix
      bool verbose{true}; ///< enable/disable verbose output
      /**
       *         Initializes all common variables.
       */
      void init();
      /**
       *         Checks if a hierarchy is suitable to the given mass spectrum.
       *         @param ho a HierarchyObject with constant isAlphab and a hierarchy candidate.
       *         @returns A bool if the hierarchy candidate is suitable to the given mass spectrum.
       */
      bool isHierarchySuitable(const HierarchyObject& ho) const;
      /**
       *         Shifts the H3m renormalization scheme to DR' scheme
       *         @param ho a HierarchyObject with constant isAlphab
       *         @return A matrix which shifts the H3m scheme to the DR' scheme at three-loop level
       *
       */
      Eigen::Matrix2d shiftH3mToDRbarPrime(const HierarchyObject& ho) const;
      /**
       *         Shifts the H3m renormalization scheme to DR' scheme
       *         @param ho a HierarchyObject with constant isAlphab
       *         @param omitLogs a flag to omit logarithmic contributions: (0) omit logs, (1) add them
       *         @return A double which shifts the H3m scheme to the DR' scheme at three-loop level
       *
       */
      double shiftH3mToDRbarPrimeMh2(const himalaya::HierarchyObject& ho, int omitLogs) const;
      std::map<unsigned int, unsigned int> flagMap{}; ///< A map which holds all hierarchy key value pairs.
      /**
       *         Maps a hierarchy to it's mother hierarchy.
       *         @param hierarchy the key to a hierarchy.
       *         @throws runtime_error Throws a runtime_error if the given hierarchy is not included.
       *         @returns The key of the mother hierarchy.
       */
      int getCorrectHierarchy(const int hierarchy) const;
      /**
       *         Fills in delta_lambda @ 3L to the given HierarchyObject
       *         @param ho a HierrachyObject
       *         @param omitXtOrders a bool to omit xtOrders of delta_lambda_EFT
       */
      void calcDeltaLambda3L(HierarchyObject& ho, bool omitXtOrders) const;
      /**
       *         Prints out some information about Himalaya.
       */
      void printInfo() const;
  };
}        // himalaya

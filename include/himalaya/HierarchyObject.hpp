// ====================================================================
// This file is part of Himalaya.
//
// Himalaya is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <array>
#include <iosfwd>
#include <string>
#include <Eigen/Core>

/**
 * @file HierarchyObject.hpp
 *
 * @brief Definition of the HierarchyObject, which contains all
 * the calculational results.
 */

namespace himalaya {
   /**
    * The HierarchyObject class.
    */
   class HierarchyObject{
   public:
      /**
       * A constructor.
       * @param isAlphab the boolean which determines wether the members are proportinal to alpha_b or alpha_t.
       */
      HierarchyObject(bool isAlphab);
      /**
       * @return the value of isAlphab.
       */
      bool getIsAlphab() const;
      /**
       * @return The key to the suitable hierarchy
       */
      int getSuitableHierarchy() const;
      /**
       * @return The absolute difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
       */
      double getAbsDiff2L() const;
      /**
       * @return The relative difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
       */
      double getRelDiff2L() const;
      /**
       * @param loops an integer which can be 1, 2 or 3.
       * @return A double which is the expansion uncertainty for the given loop order.
       */
      double getDMhExpUncertainty(int loops) const;
      /**
       * @param loops an integer which can be 0, 1, 2, 3. Here, 0 corresponds to the tree-level matrix.
       * @return The CP-even Higgs mass matrix at the given loop order.
       */
      Eigen::Matrix2d getDMh(int loops) const;
      /**
       * @param loops an integer which can be 0, 1, 2, 3. Here, 0 corresponds to the tree-level matrix.
       * @return The correction to squared CP-even Higgs mass at the given loop order.
       */
      double getDMh2(int loops) const;
      /**
       * @return The matrix M(MDR') - M(DR') at the order O(alpha_x + alpha_x*alpha_s)
       */
      Eigen::Matrix2d getDMhDRbarPrimeToMDRbarPrimeShift() const;
      /**
       * @return A vector of the MDR stop/sbottom masses. The 0th entry corresponds to the lighter particle.
       */
      Eigen::Vector2d getMDRMasses() const;
      /**
       * @return the MDRFlag integer.
       */
      int getMDRFlag() const;
      /**
       * @return Renormalization scheme key
       */
      int getRenormalizationScheme() const;
      /**
       * @return 3-loop delta_lambda with H3m logs
       */
      double getDLambdaH3m() const;
      /**
       * @return 3-loop delta_lambda with EFT logs
       */
      double getDLambdaEFT() const;
      /**
       * @return 3-loop non-logarithmic part of delta_lambda
       */
      double getDLambdaNonLog() const;
      /**
       * @return the uncertainty of delta_lambda
       * @param loops an integer, could be 0 (tree), 1 (1L), ..., 3 (3L)
       */
      double getDLambdaUncertainty(int loops) const;
      /**
       * @return delta_lambda
       * @param loops an integer, could be 0 (tree), 1 (1L), ..., 3 (3L)
       */
      double getDLambda(int loops) const;
      /**
       * @return shift to convert delta_lambda from DR' to MS scheme. This shift has to be added to the 1L result
       * @param loops an integer, could be 0 (tree), 1 (1L), ..., 3 (3L)
       */
      double getDLambdaDRbarPrimeToMSbarShift(int loops) const;
      /**
       * @return Delta_Mh2_EFT (only contributions that go with αt)
       * @param loops an integer, could be 0 (tree), 1 (1L), ..., 3 (3L)
       */
      double getDMh2EFTAt(int loops) const;
      /**
       * @return Delta_Mh2_FO
       * @param loops an integer, could be 0 (tree), 1 (1L), ..., 3 (3L)
       */
      double getDMh2FO(int loops) const;
      /**
       * @return Delta_Mh2_FO (only contributions that go with αt)
       * @param loops an integer, could be 0 (tree), 1 (1L), ..., 3 (3L)
       */
      double getDMh2FOAt(int loops) const;
      /**
       * Sets the suitable hierarchy
       * @param hierarchy the integer key of the hierarchy.
       */
      void setSuitableHierarchy(int hierarchy);
      /**
       * Sets the absolute difference of the Higgs masses at two-loop level
       * @param absDiff2L the absolute difference of the Higgs masses as a double.
       */
      void setAbsDiff2L(double absDiff2L);
      /**
       * Sets the relative difference ot the Higgs masses at two-loop level
       * @param relDiff2L the relative difference of the Higgs masses as a double.
       */
      void setRelDiff2L(double relDiff2L);
      /**
       * Sets the uncertainty of the expansion at a given loop level.
       * @param loops the integer value of the corresponding loops. Can be 1, 2 or 3.
       * @param uncertainty the expansion untertainty at the given loop order as a double.
       */
      void setDMhExpUncertainty(int loops, double uncertainty);
      /**
       * Sets the DR' -> MDR' shift
       * @param mdrShift the DR' -> MDR' shiftet matrix of the form M(MDR') - M(DR').
       */
      void setDMhDRbarPrimeToMDRbarPrimeShift(const Eigen::Matrix2d& mdrShift);
      /**
       * Sets the MDR masses
       * @param mdrMasses a vector containting the MDR masses with the lightest particle at position 0.
       */
      void setMDRMasses(const Eigen::Vector2d& mdrMasses);
      /**
       * Sets the delta of the CP-even Higgs mass matrix
       * @param loops the integer value of the corresponding loops. Can be 0, 1, 2 or 3. 0 corresponds to the tree-level.
       * @param dMh the delta of the mass matrix.
       */
      void setDMh(int loops, const Eigen::Matrix2d& dMh);
      /**
       * Sets the delta of the squared CP-even Higgs mass
       * @param loops the integer value of the corresponding loops. Can be 0, 1, 2 or 3. 0 corresponds to the tree-level.
       * @param dMh2 the delta of the squared mass.
       */
      void setDMh2(int loops, double dMh2);
      /**
       * Sets the mdrFlag to calculate the corretions in the without MDR (0) or with MDR (1) shifts
       * @param mdrFlag an int. (0) for H3m (DR')- and (1) for MDR-scheme.
       * @throws runtime_exception if the flag is neither 0 or 1 an exception is thrown.
       */
      void setMDRFlag(int mdrFlag);
      /**
       * Sets the renormalization scheme accodring to the RenScheme enum
       * @param renScheme an int according to the RenScheme enum.
       * @throws runtime_exception if the flag is not in {0,1,2,3} an exception is thrown
       */
      void setRenormalizationScheme(int renScheme);
      /**
       * Sets the Delta_lambda at loops-loop
       * @param loops an integer, could be 0 (tree), ..., 3 (3L Delta_lambda_EFT)
       * @param deltaLambda delta_lambda at tree-level
       */
      void setDLambda(int loops, double deltaLambda);
      /**
       * Sets the DR' to MS shift for delta_lambda at loops-loop
       * @param loops an integer, could be 0(tree), ..., 3 (3L Delta_lambda_EFT shift)
       * @param shift the shift
       */
      void setDLambdaDRbarPrimeToMSbarShift(int loops, double shift);
      /**
       * Sets the delta_lambda at 3-loop order with H3m logs.
       * This variable is only used in the hierarchy selection process.
       * @param deltaLambda delta_lambda at 3-loop order.
       */
      void setDLambdaH3m(double deltaLambda);
      /**
       * Sets the delta_lambda at 3-loop order with EFT logs
       * @param deltaLambda delta_lambda at 3-loop order.
       */
      void setDLambdaEFT(double deltaLambda);
      /**
       * Sets the non-logarithmic part of delta_lambda at 3-loop order.
       * @param deltaLambda constant part of delta_lambda at 3-loop order.
       */
      void setDLambdaNonLog(double deltaLambda);
      /**
       * Sets the Xt parts of the uncertainty of delta_lambda_EFT
       * @param uncertainty of 3-loop delta_lambda
       */
      void setDLambdaXtUncertainty(double uncertainty);
      /**
       * Returns the H3m notation of a given hierarchy.
       * @param hierarchy An integer of a Himalaya hierarchy.
       * @return Returns the corresponding H3m notation of the given hierarchy as a string.
       */
      std::string getH3mHierarchyNotation(int hierarchy) const;
      /**
       * Sets the DR' -> H3m shift which should be added to the DR' result
       * @param shift the DR' -> H3m shift
       */
      void setDMhDRbarPrimeToH3mShift(const Eigen::Matrix2d& shift);
      /**
       * Returns the DR' -> H3m shift which should be added to the DR' result
       * @return Returns the DR' -> H3m shift
       */
      Eigen::Matrix2d getDMhDRbarPrimeToH3mShift() const;
      /**
       * Sets Delta_Mh2_EFT at loops-loop (only contributions that go with αt)
       * @param loops an integer, could be 0 (tree), ..., 3 (3L with Delta_lambda_EFT)
       * @param deltaMh2 delta_Mh^2
       */
      void setDMh2EFT(int loops, double deltaMh2);
      /**
       * Sets Delta_Mh2_FO at loops-loop
       * @param loops an integer, could be 0 (tree), ..., 3 (3L)
       * @param deltaMh2 delta_Mh^2
       */
       void setDMh2FO(int loops, double deltaMh2);
      /**
       * Sets Delta_Mh2_FO at loops-loop (only contributions that go with αt)
       * @param loops an integer, could be 0 (tree), ..., 3 (3L)
       * @param deltaMh2 delta_Mh^2
       */
       void setDMh2FOAt(int loops, double deltaMh2);
   private:
      bool isAlphab{false};                    ///< the bool isAlphab
      int hierarchy{};                         ///< the suitable hierarchy
      int mdrFlag{};                           ///< the MDR-scheme flag
      int renormalizationScheme{};             ///< the renormalization scheme flag
      double absDiff2L{};                      ///< the absolute difference of the two loop Higgs masses
      double relDiff2L{};                      ///< the relative difference of the two loop Higgs masses
      std::array<double,4> expUncertainties{}; ///< holds the expansion uncertainties, the keys are the loop order: 0, 1, 2, 3
      std::array<Eigen::Matrix2d,4> dMhMap{ {Eigen::Matrix2d::Zero(),
                                             Eigen::Matrix2d::Zero(),
                                             Eigen::Matrix2d::Zero(),
                                             Eigen::Matrix2d::Zero()} }; ///< holds all mass matrices at the given loop order
      std::array<double,4> dMh2Map{};          ///< holds loop corretions to squared CP-even Higgs mass
      Eigen::Matrix2d mdrShift{Eigen::Matrix2d::Zero()};  ///< the mass matrix of the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s
      Eigen::Vector2d mdrMasses{Eigen::Vector2d::Zero()}; ///< the 'array' which holds the MDR masses
      double dLambdaH3m{};                     ///< delta_lambda 3-loop with H3m logs
      double dLambdaEFT{};                     ///< delta_lambda 3-loop with EFT logs
      double dLambdaNonLog{};                  ///< delta_lambda 3-loop, non-logarithmic part
      double dLambdaXtUncertainty{};           ///< The uncertainty of delta_lambda_EFT due to mising Xt terms
      Eigen::Matrix2d h3mShift{Eigen::Matrix2d::Zero()}; ///< The DR' -> H3m shift which should be added to the DR' result
      std::array<double,4> dLambdaMap{};       ///< holds all delta_lambda corrections multiplied with prefactors
      std::array<double,4> dLambdaDRbarPrimeToMSbarShiftMap{}; ///< holds all DR' -> MS shifts for delta_lambda corrections multiplied with prefactors
      std::array<double,4> dMh2EFTMap{};       ///< holds all delta_Mh2_EFT corrections
      std::array<double,4> dMh2FOMap{};        ///< holds all delta_Mh2_FO corrections
      std::array<double,4> dMh2FOAtMap{};      ///< holds all delta_Mh2_FO corrections (only contributions that go with αt)
   };
   /**
    * Prints out all information of the HierarchyObject
    */
   std::ostream& operator<<(std::ostream&, const HierarchyObject&);
}        // himalaya

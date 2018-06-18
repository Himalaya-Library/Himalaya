#pragma once

#include "version.hpp"
#include "Hierarchies.hpp"
#include <iosfwd>
#include <Eigen/Eigenvalues>
#include <vector>
#include <map>

namespace himalaya{
   /**
    * 	The HierarchyObject class.
    */
   class HierarchyObject{
   public:
      /**
       * 	A constructor.
       * 	@param isAlphab the boolean which determines wether the members are proportinal to alpha_b or alpha_t.
       */
      HierarchyObject(bool isAlphab);
      /**
       * 	@return the value of isAlphab.
       */
      bool getIsAlphab() const;
      /**
       * 	@return The key to the suitable hierarchy
       */
      int getSuitableHierarchy() const;
      /**
       * 	@return The absolute difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
       */
      double getAbsDiff2L() const;
      /**
       * 	@return The relative difference of the exact and expanded Higgs masses at two-loop level at the order O(alpha_x + alpha_x*alpha_s).
       */
      double getRelDiff2L() const;
      /**
       * 	@param loops an integer which can be 1, 2 or 3.
       * 	@return A double which is the expansion uncertainty for the given loop order.
       */
      double getExpUncertainty(int loops) const;
      /**
       * 	@param loops an integer which can be 0, 1, 2, 3. Here, 0 corresponds to the tree-level matrix.
       * 	@return The CP-even Higgs mass matrix at the given loop order.
       */
      Eigen::Matrix2d getDMh(int loops) const;
      /**
       * 	@return The matrix M(MDR') - M(DR') at the order O(alpha_x + alpha_x*alpha_s)
       */
      Eigen::Matrix2d getDRbarPrimeToMDRbarPrimeShift() const;
      /**
       * 	@return A vector of the MDR stop/sbottom masses. The 0th entry corresponds to the lighter particle.
       */
      Eigen::Vector2d getMDRMasses() const;
      /**
       * 	@return the MDRFlag integer.
       */
      int getMDRFlag() const;
      /**
       * 	@return Renormalization scheme key
       */
      int getRenormalizationScheme() const;
      /**
       * 	@return 3-loop delta_lambda with H3m logs
       */
      double getDeltaLambdaH3m() const;
      /**
       * 	@return 3-loop delta_lambda with EFT logs
       */
      double getDeltaLambdaEFT() const;
      /**
       *        @return 3-loop non-logarithmic part of delta_lambda
       */
      double getDeltaLambdaNonLog() const;
      /**
       *        @return uncertainty of 3-loop delta_lambda_H3m
       */
      double getDeltaLambdaUncertaintyH3m() const;
      /**
       *        @return uncertainty of 3-loop delta_lambda_EFT
       */
      double getDeltaLambdaUncertaintyEFT() const;
      /**
       * 	@return the DR' -> MS shift for delta_lambda_H3m which should be added to the DR' result
       */
      double getDRbarPrimeToMSbarShiftH3m() const;
      /**
       * 	@return the DR' -> MS shift for delta_lambda_EFT which should be added to the DR' result
       */
      double getDRbarPrimeToMSbarShiftEFT() const;
      /**
       * 	@return delta_lambda @ tree-level
       */
      double getDeltaLambda0L() const;
      /**
       * 	@return delta_lambda @ 1-loop order
       */
      double getDeltaLambda1L() const;
      /**
       * 	@return delta_lambda @ 2-loop order
       */
      double getDeltaLambda2L() const;
      /**
       * 	@return shift to convert 1L delta_lambda from DR' to MS scheme. This shift has to be added to the 1L result
       */
      double getDRbarPrimeToMSbarShiftDeltaLambda1L() const;
      /**
       * 	@return shift to voncert 2L delta_lambda from DR' to MS scheme. This shift has to be added to the 2L result
       */
      double getDRbarPrimeToMSbarShiftDeltaLambda2L() const;
      /**
       * 	@return Delta_Mh2_EFT @ 1L
       */
      double getDeltaMh2EFT1L() const;
      /**
       * 	@return Delta_Mh2_EFT @ 2L
       */
      double getDeltaMh2EFT2L() const;
      /**
       * 	@return Delta_Mh2_EFT @ 3L
       */
      double getDeltaMh2EFT3L() const;
      /**
       * 	Sets the suitable hierarchy
       * 	@param hierarchy the integer key of the hierarchy.
       */
      void setSuitableHierarchy(int hierarchy);
      /**
       * 	Sets the absolute difference of the Higgs masses at two-loop level
       * 	@param absDiff2L the absolute difference of the Higgs masses as a double.
       */
      void setAbsDiff2L(double absDiff2L);
      /**
       * 	Sets the relative difference ot the Higgs masses at two-loop level
       * 	@param relDiff2L the relative difference of the Higgs masses as a double.
       */
      void setRelDiff2L(double relDiff2L);
      /**
       * 	Sets the uncertainty of the expansion at a given loop level.
       * 	@param loops the integer value of the corresponding loops. Can be 1, 2 or 3.
       * 	@param uncertainty the expansion untertainty at the given loop order as a double.
       */
      void setExpUncertainty(int loops, double uncertainty);
      /**
       * 	Sets the DR' -> MDR' shift
       * 	@param mdrShift the DR' -> MDR' shiftet matrix of the form M(MDR') - M(DR').
       */
      void setDRbarPrimeToMDRbarPrimeShift(const Eigen::Matrix2d& mdrShift);
      /**
       * 	Sets the MDR masses
       * 	@param mdrMasses a vector containting the MDR masses with the lightest particle at position 0.
       */
      void setMDRMasses(Eigen::Vector2d& mdrMasses);
      /**
       * 	Sets the delta of the CP-even Higgs mass matrix
       * 	@param loops the integer value of the corresponding loops. Can be 0, 1, 2 or 3. 0 corresponds to the tree-level.
       * 	@param dMh the delta of the mass matrix.
       */
      void setDMh(int loops, const Eigen::Matrix2d& dMh);
      /**
       * 	Sets the mdrFlag to calculate the corretions in the without MDR (0) or with MDR (1) shifts
       * 	@param mdrFlag an int. (0) for H3m (DR')- and (1) for MDR-scheme.
       * 	@throws runtime_exception if the flag is neither 0 or 1 an exception is thrown.
       */
      void setMDRFlag(int mdrFlag);
      /**
       * 	Sets the renormalization scheme accodring to the RenScheme enum
       * 	@param renScheme an int according to the RenScheme enum.
       * 	@throws runtime_exception if the flag is not in {0,1,2,3} an exception is thrown
       */
      void setRenormalizationScheme(int renScheme);
      /**
       * 	Sets the delta_lambda at tree-level
       * 	@param deltaLambda delta_lambda at tree-level
       */
      void setDeltaLambda0L(double deltaLambda);
      /**
       * 	Sets the delta_lambda at 1-loop order
       * 	@param deltaLambda delta_lambda at 1-loop order.
       */
      void setDeltaLambda1L(double deltaLambda);
      /**
       * 	Sets the delta_lambda at 2-loop order
       * 	@param deltaLambda delta_lambda at 2-loop order
       */
      void setDeltaLambda2L(double deltaLambda);
      /**
       * 	Sets the DR' to MS shift for 1L delta_lambda
       * 	@param shift the shift
       */
      void setDRbarPrimeToMSbarShiftDeltaLambda1L(double shift);
      /**
       * 	Sets the DR' to MS shift for 2L delta_lambda
       * 	@param shift the shift
       */
      void setDRbarPrimeToMSbarShiftDeltaLambda2L(double shift);
      /**
       * 	Sets the delta_lambda at 3-loop order with H3m logs.
       * 	@param deltaLambda delta_lambda at 3-loop order.
       */
      void setDeltaLambdaH3m(double deltaLambda);
      /**
       * 	Sets the delta_lambda at 3-loop order with EFT logs
       * 	@param deltaLambda delta_lambda at 3-loop order.
       */
      void setDeltaLambdaEFT(double deltaLambda);
      /**
       * 	Sets the non-logarithmic part of delta_lambda at 3-loop order.
       * 	@param deltaLambda constant part of delta_lambda at 3-loop order.
       */
      void setDeltaLambdaNonLog(double deltaLambda);
      /**
       * 	Sets the Xt parts of the uncertainty of delta_lambda_H3m
       *        @param uncertainty of 3-loop delta_lambda
       */
      void setDeltaLambdaXtUncertaintyH3m(double uncertainty);
      /**
       * 	Sets the Xt parts of the uncertainty of delta_lambda_EFT
       *        @param uncertainty of 3-loop delta_lambda
       */
      void setDeltaLambdaXtUncertaintyEFT(double uncertainty);
      /**
       * 	Sets the DR' -> MS shift for delta_lambda_H3m which should be added to the DR' result
       * 	@param shift the DR' -> MS shift which should be added to the 3-loop threshold correction
       */
      void setDRbarPrimeToMSbarShiftH3m(double shift);
      /**
       * 	Sets the DR' -> MS shift for delta_lambda_EFT which should be added to the DR' result
       * 	@param shift the DR' -> MS shift which should be added to the 3-loop threshold correction
       */
      void setDRbarPrimeToMSbarShiftEFT(double shift);
      /**
       *	Returns the H3m notation of a given hierarchy.
       *	@param hierarchy An integer of a Himalaya hierarchy.
       *	@return Returns the corresponding H3m notation of the given hierarchy as a string.
       */
      std::string getH3mHierarchyNotation(int hierarchy) const;
      /**
       * 	Sets the DR' -> H3m shift which should be added to the DR' result
       * 	@param shift the DR' -> H3m shift
       */
      void setDRbarPrimeToH3mShift(const Eigen::Matrix2d shift);
      /**
       * 	Returns the DR' -> H3m shift which should be added to the DR' result
       * 	@return Returns the DR' -> H3m shift
       */
      Eigen::Matrix2d getDRbarPrimeToH3mShift() const;
      /**
       * 	Sets Delta_Mh2_EFT @ 1L O(at)
       * 	@param deltaMh2 delta_Mh^2
       */
      void setDeltaMh2EFT1L(double deltaMh2);
      /**
       * 	Sets Delta_Mh2_EFT @ 2L O(at*as)
       * 	@param deltaMh2 delta_Mh^2
       */
      void setDeltaMh2EFT2L(double deltaMh2);
      /**
       * 	Sets Delta_Mh2_EFT @ 3L O(at*as^2)
       * 	@param deltaMh2 delta_Mh^2
       */
      void setDeltaMh2EFT3L(double deltaMh2);
      /**
       * 	Set the expasion uncertainty for delta_lambda
       * 	@param expUncertLambda the expansion uncertainty for delta_lambda
       */
      void setExpUncertaintyDeltaLambda(double expUncertLambda);
   private:
      bool isAlphab{false};								/**< the bool isAlphab */
      int hierarchy{};									/**< the suitable hierarchy */
      int mdrFlag{0};									/**< the MDR-scheme flag */
      int renormalizationScheme{RenSchemes::DRBARPRIME};				/**< the renormalization scheme flag */
      double absDiff2L{};								/**< the absolute difference of the two loop Higgs masses */
      double relDiff2L{};								/**< the relative difference of the two loop Higgs masses */
      std::map<int, double> expUncertainties{};						/**< the map which holds the expansion uncertainties, the keys are the loop order: 1, 2, 3 */
      std::map<int, Eigen::Matrix2d> dMhMap{};						/**< the map which holds all mass matrices at the given loop order */
      Eigen::Matrix2d mdrShift{};							/**< the mass matrix of the difference of the MDR - DR contributions of the order alpha_x + alpha_x*alpha_s */
      Eigen::Vector2d mdrMasses{};							/**< the 'vector' which holds the MDR masses */
      double deltaLambdaH3m{};								/**< delta_lambda 3-loop with H3m logs */
      double deltaLambdaEFT{};								/**< delta_lambda 3-loop with EFT logs */
      double deltaLambdaNonLog{};                                                       /**< delta_lambda 3-loop, non-logarithmic part */
      double drBarPrimeToMSbarShiftH3m{};						/**< The shift to convert the DR' delta_lambda_H3m to the MS scheme */
      double drBarPrimeToMSbarShiftEFT{};						/**< The shift to convert the DR' delta_lambda_EFT to the MS scheme */
      double deltaLambdaXtUncertaintyH3m{};						/**< The uncertainty of delta_lambda_H3m due to mising Xt terms */
      double deltaLambdaXtUncertaintyEFT{};						/**< The uncertainty of delta_lambda_EFT due to mising Xt terms */
      double expansionUncertaintyDeltaLambda{};						/**< The expansion uncertainty of delta_lambda */
      Eigen::Matrix2d h3mShift{};							/**< The DR' -> H3m shift which should be added to the DR' result */
      double deltaLambda0L{};								/**< Delta_lambda @ 0L */
      double deltaLambda1L{};								/**< Delta_lambda @ 1L O(at) */
      double deltaLambda2L{};								/**< Delta_lambda @ 2L O(at*as) */
      double drToMSDL1L{};								/**< The shift to convert delta_lambda from DR' to MS @ 1L */
      double drToMSDL2L{};								/**< The shift to convert delta_lambda from DR' to MS @ 2L */
      double mh2EFT1L{};								/**< Mh2_EFT @ 1L O(at) */
      double mh2EFT2L{};								/**< Mh2_EFT @ 2L O(at*as) */
      double mh2EFT3L{};								/**< Mh2_EFT @ 3L O(at*as^2) */
      /**
       * 	Sorts a vector.
       * 	@param vector The vector which should be sorted.
       *	@return Returns a vector the lightest entry at position 0.
       */
      Eigen::Vector2d sortVector(Eigen::Vector2d& vector);
   };
   /**
    * 	Prints out all information of the HierarchyObject
    */
   std::ostream& operator<<(std::ostream&, const HierarchyObject&);
}	// himalaya
